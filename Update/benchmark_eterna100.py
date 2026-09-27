"""
benchmark_eterna100.py

Benchmark the information-field solver on the full Eterna100 dataset
using the ViennaRNA Turner nearest-neighbor model.

Success criteria:
  - MFE:  MFE(sequence) == target_structure  (strict)
  - uMFE: MFE(sequence) == target_structure AND target is unique MFE
  - Relaxed: base-pair distance <= 3  (near-correct)

Outputs:
  results/benchmark_eterna100_<timestamp>.log   full transcript
  results/benchmark_eterna100.csv               per-puzzle data

Layered statistics are reported by sequence length:
  Short   : L <= 100
  Medium  : 100 < L <= 200
  Long    : L > 200

Usage:
  python benchmark_eterna100.py
"""

import os
import sys
import csv
import time
from datetime import datetime
from multiprocessing import Pool

# Ensure single-threaded ViennaRNA per worker
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import RNA  # loud failure if ViennaRNA not installed
from rna_common import (
    metropolis_sample, structure_distance, fold_mfe,
    parse_structure_features, load_eterna100_from_url,
)


# ============================================================
# Layered statistics
# ============================================================

def length_bin(L):
    if L <= 100:
        return "short (L<=100)"
    elif L <= 200:
        return "medium (100<L<=200)"
    else:
        return "long (L>200)"


def helix_bin(H):
    if H <= 3:
        return "simple (H<=3)"
    elif H <= 7:
        return "moderate (4<=H<=7)"
    else:
        return "complex (H>=8)"


# ============================================================
# Worker
# ============================================================

def _run_puzzle(args):
    """Run MH sampling on one puzzle; return per-puzzle result."""
    name, structure, T, seed, n_steps, n_samples, burn_in = args

    t0 = time.time()
    r = metropolis_sample(
        structure, T=T, lambda_dist=5.0,
        n_steps=n_steps, n_samples=n_samples, burn_in=burn_in,
        seed=seed,
    )
    elapsed = time.time() - t0

    # Evaluate all samples; keep the best (closest to target)
    samples = r["samples"]
    best_dist = 999
    best_seq = None
    best_mfe = None
    for s in samples:
        mfe, _ = fold_mfe(s)
        d = structure_distance(mfe, structure)
        if d < best_dist:
            best_dist = d
            best_seq = s
            best_mfe = mfe
            if d == 0:
                break

    # Criteria
    mfe_solved = (best_dist == 0)
    relaxed_solved = (best_dist <= 3)

    # uMFE: target is the unique MFE
    u_mfe_solved = False
    if mfe_solved and best_seq is not None:
        try:
            fc = RNA.fold_compound(best_seq, None, RNA.OPTION_EVAL_ONLY)
            target_energy = fc.eval_structure(structure)
            # subopt within 1 kcal/mol
            subopt = RNA.subopt(best_seq, 1.0)
            n_same = sum(1 for _, e in subopt
                         if abs(e - target_energy) < 1e-6)
            u_mfe_solved = (n_same == 1)
        except Exception:
            u_mfe_solved = False

    # Structural features
    feats = parse_structure_features(structure)

    return {
        "name": name,
        "length": feats["length"],
        "num_helices": feats["num_helices"],
        "max_loop": feats["max_loop"],
        "mfe_solved": int(mfe_solved),
        "u_mfe_solved": int(u_mfe_solved),
        "relaxed_solved": int(relaxed_solved),
        "best_dist": best_dist,
        "elapsed_s": round(elapsed, 2),
        "length_bin": length_bin(feats["length"]),
        "helix_bin": helix_bin(feats["num_helices"]),
    }


# ============================================================
# Logging helper
# ============================================================

class Tee:
    def __init__(self, path):
        self.terminal = sys.stdout
        self.logfile = open(path, "w", encoding="utf-8")

    def write(self, msg):
        self.terminal.write(msg)
        self.logfile.write(msg)
        self.logfile.flush()

    def flush(self):
        self.terminal.flush()
        self.logfile.flush()


# ============================================================
# Statistics printer
# ============================================================

def print_layered(rows, key, write):
    """Print success rates for each bin of `key`."""
    bins = {}
    for r in rows:
        b = r[key]
        if b not in bins:
            bins[b] = []
        bins[b].append(r)

    write(f"\nLayered by {key}:")
    write(f"  {'bin':<26}{'n':>4}{'MFE':>10}{'uMFE':>10}"
          f"{'Relaxed':>10}{'MeanDist':>10}")
    write("  " + "-" * 68)
    for b in sorted(bins.keys()):
        sub = bins[b]
        n = len(sub)
        mfe = sum(r["mfe_solved"] for r in sub)
        umfe = sum(r["u_mfe_solved"] for r in sub)
        rel = sum(r["relaxed_solved"] for r in sub)
        md = sum(r["best_dist"] for r in sub) / n
        write(f"  {b:<26}{n:>4}"
              f"{mfe:>4}/{n:<3} ({100*mfe/n:>5.1f}%)"
              f"  "
              f"{umfe:>3}/{n:<3} ({100*umfe/n:>5.1f}%)"
              f"  "
              f"{rel:>3}/{n:<3} ({100*rel/n:>5.1f}%)"
              f"  {md:>8.2f}")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    log_path = f"results/benchmark_eterna100_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    sys.stdout = Tee(log_path)

    write = sys.stdout.write

    # ---------- Header ----------
    write("=" * 70 + "\n")
    write("Eterna100 Benchmark (Full Turner Model)\n")
    write("=" * 70 + "\n")
    write(f"Timestamp  : {datetime.now().isoformat(timespec='seconds')}\n")
    write(f"ViennaRNA  : {RNA.__version__}\n")

    # ---------- Parameters ----------
    T = 1.0
    SEED = 42
    N_STEPS = 3000
    N_SAMPLES = 300
    BURN_IN = 800
    N_WORKERS = 2

    write(f"Temperature: {T}\n")
    write(f"Seed       : {SEED}\n")
    write(f"MH steps   : {N_STEPS} (burn-in {BURN_IN}, samples {N_SAMPLES})\n")
    write(f"Workers    : {N_WORKERS}\n")
    write(f"Lambda     : 5.0\n\n")

    # ---------- Load puzzles (NO length cap) ----------
    puzzles = load_eterna100_from_url()
    puzzles.sort(key=lambda p: len(p["structure"]))
    selected = puzzles  # all valid puzzles

    write(f"Loaded {len(selected)} puzzles (no length cap)\n")
    if selected:
        write(f"Length range: {len(selected[0]['structure'])} – "
              f"{len(selected[-1]['structure'])} nt\n\n")

    write("Full puzzle list:\n")
    for i, p in enumerate(selected, 1):
        feats = parse_structure_features(p["structure"])
        write(f"  [{i:>2}] L={feats['length']:>4}  "
              f"H={feats['num_helices']:>2}  "
              f"{p['title'][:55]}\n")
    write("\n")

    # ---------- Run ----------
    write("=" * 70 + "\n")
    write("Running benchmark\n")
    write("=" * 70 + "\n")

    start = time.time()
    rows = []
    args = [(p["title"], p["structure"], T, SEED,
             N_STEPS, N_SAMPLES, BURN_IN) for p in selected]

    with Pool(processes=N_WORKERS) as pool:
        for i, r in enumerate(pool.imap(_run_puzzle, args), 1):
            rows.append(r)
            write(f"[{i:>2}/{len(selected)}] "
                  f"{r['name'][:32]:<32} "
                  f"L={r['length']:>4}  H={r['num_helices']:>2}  "
                  f"MFE={'OK' if r['mfe_solved'] else '-':>2}  "
                  f"uMFE={'OK' if r['u_mfe_solved'] else '-':>2}  "
                  f"rel={'OK' if r['relaxed_solved'] else '-':>2}  "
                  f"dist={r['best_dist']:>2}  "
                  f"({r['elapsed_s']}s)\n")

    total = time.time() - start

    # ---------- Summary ----------
    n = len(rows)
    n_mfe = sum(r["mfe_solved"] for r in rows)
    n_u_mfe = sum(r["u_mfe_solved"] for r in rows)
    n_rel = sum(r["relaxed_solved"] for r in rows)
    mean_dist = sum(r["best_dist"] for r in rows) / n if n else 0

    write("\n" + "=" * 70 + "\n")
    write("Overall summary\n")
    write("=" * 70 + "\n")
    write(f"Puzzles tested     : {n}\n")
    write(f"MFE solved         : {n_mfe}/{n} ({100*n_mfe/n:.1f}%)\n")
    write(f"uMFE solved        : {n_u_mfe}/{n} ({100*n_u_mfe/n:.1f}%)\n")
    write(f"Relaxed (dist<=3)  : {n_rel}/{n} ({100*n_rel/n:.1f}%)\n")
    write(f"Mean best distance : {mean_dist:.2f}\n")
    write(f"Total wall time    : {total:.1f} s ({total/60:.1f} min)\n")
    write(f"Mean per puzzle    : {total/n:.2f} s\n")

    # ---------- Layered statistics ----------
    write("\n" + "=" * 70 + "\n")
    write("Layered statistics\n")
    write("=" * 70)

    print_layered(rows, "length_bin", write)
    print_layered(rows, "helix_bin", write)

    # ---------- Save CSV ----------
    csv_path = "results/benchmark_eterna100.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    write(f"\n\nWrote {csv_path}\n")
    write(f"Wrote {log_path}\n")

    # ---------- Comparison with published methods ----------
    write("\n" + "=" * 70 + "\n")
    write("Comparison with published methods (MFE criterion)\n")
    write("=" * 70 + "\n")
    write(f"{'Method':<32}{'MFE Solved':>16}{'Runtime':>14}\n")
    write("-" * 62 + "\n")
    write(f"{'DesiRNA (2025)':<32}{'100/100 (100%)':>16}{'24 h':>14}\n")
    write(f"{'This work (MFE, strict)':<32}"
          f"{f'{n_mfe}/{n} ({100*n_mfe/n:.0f}%)':>16}"
          f"{f'{total/60:.1f} min':>14}\n")
    write(f"{'This work (relaxed, dist<=3)':<32}"
          f"{f'{n_rel}/{n} ({100*n_rel/n:.0f}%)':>16}"
          f"{f'{total/60:.1f} min':>14}\n")

    # ---------- Close ----------
    sys.stdout = sys.stdout.terminal
    if hasattr(sys.stdout, "logfile"):
        pass