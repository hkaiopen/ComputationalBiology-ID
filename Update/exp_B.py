"""
Experiment B: Structural dependence of the correct-folding basin.

For 15 short puzzles (L <= 60), run 100 seeds at T=1.0.
Correlate basin fraction with helix count, max loop, and length.
"""

import sys
import csv
import time
from multiprocessing import Pool

from rna_common import (
    metropolis_sample, parse_structure_features,
    load_eterna100_from_url, make_log, RNA,
)


def _worker(args):
    structure, T, seed = args
    r = metropolis_sample(structure, T=T, lambda_dist=5.0,
                          n_steps=3000, n_samples=300, burn_in=800,
                          seed=seed)
    return seed, r["exact_match"], r["avg_distance"]


def parallel_seeds(structure, T, seeds, n_workers):
    with Pool(processes=n_workers) as pool:
        return pool.map(_worker, [(structure, T, s) for s in seeds])


if __name__ == "__main__":
    log_path, tee = make_log("exp_B")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    puzzles.sort(key=lambda p: len(p["structure"]))
    selected = [p for p in puzzles if len(p["structure"]) <= 60][:15]

    print("=" * 70)
    print("EXPERIMENT B: Structural dependence")
    print("=" * 70)
    print(f"Selected {len(selected)} puzzles (L <= 60)\n")

    SEEDS = list(range(42, 142))
    N_WORKERS = 2
    rows = []
    start = time.time()

    for idx, p in enumerate(selected, 1):
        feats = parse_structure_features(p["structure"])
        t0 = time.time()
        results = parallel_seeds(p["structure"], 1.0, SEEDS, N_WORKERS)
        exacts = [r[1] for r in results]
        dists = [r[2] for r in results]
        basin = sum(1 for d in dists if d < 1.0) / len(dists)

        row = {
            "name": p["title"][:40],
            "length": feats["length"],
            "num_helices": feats["num_helices"],
            "max_loop": feats["max_loop"],
            "paired_ratio": round(feats["paired_ratio"], 3),
            "mean_exact": round(sum(exacts) / len(exacts), 4),
            "basin_fraction": round(basin, 4),
            "mean_dist": round(sum(dists) / len(dists), 2),
        }
        rows.append(row)
        print(f"[{idx:>2}/{len(selected)}] {p['title'][:28]:<28} "
              f"L={feats['length']:>3} H={feats['num_helices']:>2} "
              f"exact={row['mean_exact']:>6.1%} "
              f"basin={row['basin_fraction']:>6.1%} "
              f"({time.time()-t0:.0f}s)")

        with open("results/exp_B_structural_dependence.csv",
                  "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    print(f"\nTotal: {time.time()-start:.0f}s")
    sys.stdout = tee.terminal
    tee.logfile.close()