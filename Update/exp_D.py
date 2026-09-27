"""
Experiment D: Sampler robustness.

Compare three MH variants:
  (1) single-site mutation
  (2) block mutation (mutate a helix as a unit)
  (3) simulated annealing (temperature cools over time)

100 seeds each, T=1.0.
"""

import sys
import math
import random
import csv
import time

from rna_common import (
    dual_energy, structure_distance, fold_mfe,
    _find_partner, _extract_pairs,
    load_eterna100_from_url, make_log, RNA,
)


def propose_block(seq, structure, rng):
    """Mutate a random helix (contiguous paired segment) as a unit."""
    pairs = _extract_pairs(structure)
    if not pairs:
        return seq, False
    i, j = rng.choice(list(pairs))
    n = len(seq)
    # Expand to contiguous paired positions around i
    start = i
    while start - 1 >= 0 and structure[start - 1] in '()':
        start -= 1
    end = i
    while end + 1 < n and structure[end + 1] in '()':
        end += 1

    seq_list = list(seq)
    changed = False
    for pos in range(start, end + 1):
        if structure[pos] in '()':
            partner = _find_partner(structure, pos)
            if partner is None:
                continue
            pb = seq_list[partner]
            if pb == 'A':
                cands = ['U']
            elif pb == 'U':
                cands = ['A', 'G']
            elif pb == 'G':
                cands = ['C', 'U']
            elif pb == 'C':
                cands = ['G']
            else:
                cands = ['A', 'U', 'G', 'C']
            new = rng.choice(cands)
            if new != seq_list[pos]:
                seq_list[pos] = new
                changed = True
    return ''.join(seq_list), changed


def propose_single(seq, structure, rng):
    from rna_common import propose_mutation
    return propose_mutation(seq, structure, rng)


def run_sampler(structure, T, lambda_dist, seed, mode, n_steps=3000,
                n_samples=300, burn_in=800):
    rng = random.Random(seed)
    n = len(structure)
    seq = ''.join(rng.choice('AUGC') for _ in range(n))
    curr_E = dual_energy(seq, structure, lambda_dist)
    samples = []
    interval = max(1, (n_steps - burn_in) // n_samples)

    for step in range(n_steps):
        if mode == "annealing":
            T_eff = T * (2.0 - step / n_steps)
        else:
            T_eff = T

        if mode == "block":
            new_seq, changed = propose_block(seq, structure, rng)
        else:
            new_seq, changed = propose_single(seq, structure, rng)

        if changed:
            new_E = dual_energy(new_seq, structure, lambda_dist)
            dE = new_E - curr_E
            if dE <= 0 or rng.random() < math.exp(-dE / T_eff):
                seq = new_seq
                curr_E = new_E
        if step >= burn_in and (step - burn_in) % interval == 0:
            samples.append(seq)

    exact = 0
    dists = []
    for s in samples:
        mfe, _ = fold_mfe(s)
        d = structure_distance(mfe, structure)
        dists.append(d)
        if d == 0:
            exact += 1
    return {
        "exact_match": exact / len(samples) if samples else 0.0,
        "avg_distance": sum(dists) / len(dists) if dists else 0.0,
    }


if __name__ == "__main__":
    log_path, tee = make_log("exp_D")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    target = next(p for p in puzzles
                  if p["title"].strip().startswith("[CloudBeta] An Arm and a Leg"))
    structure = target["structure"]

    print("=" * 70)
    print("EXPERIMENT D: Sampler robustness")
    print("=" * 70)
    print(f"Puzzle: {target['title']}\n")

    SEEDS = list(range(42, 142))
    modes = ["single-site", "block", "annealing"]

    rows = []
    for mode in modes:
        print(f"\nSampler: {mode}")
        t0 = time.time()
        exacts, dists = [], []
        for seed in SEEDS:
            r = run_sampler(structure, 1.0, 5.0, seed, mode)
            exacts.append(r["exact_match"])
            dists.append(r["avg_distance"])
        basin = sum(1 for d in dists if d < 1.0) / len(dists)
        row = {
            "sampler": mode,
            "mean_exact": round(sum(exacts) / len(exacts), 4),
            "basin_fraction": round(basin, 4),
            "mean_dist": round(sum(dists) / len(dists), 2),
        }
        rows.append(row)
        print(f"  mean_exact={row['mean_exact']:.2%}  "
              f"basin={row['basin_fraction']:.1%}  "
              f"dist={row['mean_dist']:.2f}  ({time.time()-t0:.0f}s)")

    with open("results/exp_D_sampler_robustness.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote results/exp_D_sampler_robustness.csv")
    sys.stdout = tee.terminal
    tee.logfile.close()