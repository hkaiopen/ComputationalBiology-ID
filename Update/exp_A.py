"""
Experiment A: Temperature stability of the correct-folding basin.

For a single puzzle, run 100 seeds at each of 5 temperatures.
Report correct-basin fraction (dist < 1) and mean structure distance.
"""

import os
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
    log_path, tee = make_log("exp_A")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    target = next(p for p in puzzles
                  if p["title"].strip().startswith("[CloudBeta] An Arm and a Leg"))
    structure = target["structure"]
    feats = parse_structure_features(structure)

    print("=" * 70)
    print("EXPERIMENT A: Temperature stability")
    print("=" * 70)
    print(f"Puzzle: {target['title']}")
    print(f"Length: {feats['length']}, helices: {feats['num_helices']}\n")

    TEMPERATURES = [0.5, 0.75, 1.0, 1.5, 2.0]
    SEEDS = list(range(42, 142))
    N_WORKERS = 2

    rows = []
    per_seed_T1 = None
    start = time.time()

    for T in TEMPERATURES:
        t0 = time.time()
        results = parallel_seeds(structure, T, SEEDS, N_WORKERS)
        exacts = [r[1] for r in results]
        dists = [r[2] for r in results]
        basin = sum(1 for d in dists if d < 1.0) / len(dists)

        row = {
            "T": T,
            "mean_exact": round(sum(exacts) / len(exacts), 4),
            "basin_fraction": round(basin, 4),
            "mean_dist": round(sum(dists) / len(dists), 2),
        }
        rows.append(row)
        if T == 1.0:
            per_seed_T1 = results

        print(f"T={T:.2f}  exact={row['mean_exact']:.2%}  "
              f"basin={row['basin_fraction']:.1%}  "
              f"dist={row['mean_dist']:.2f}  ({time.time()-t0:.0f}s)")

    print(f"\nTotal: {time.time()-start:.0f}s")

    with open("results/exp_A_temperature_stability.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["T", "mean_exact",
                                          "basin_fraction", "mean_dist"])
        w.writeheader()
        w.writerows(rows)
    print("Wrote results/exp_A_temperature_stability.csv")

    if per_seed_T1:
        with open("results/exp_A_T1_per_seed.csv",
                  "w", newline="", encoding="utf-8") as f:
            f.write("seed,exact,avg_dist\n")
            for seed, exact, dist in per_seed_T1:
                f.write(f"{seed},{exact:.4f},{dist:.2f}\n")
        print("Wrote results/exp_A_T1_per_seed.csv")

    sys.stdout = tee.terminal
    tee.logfile.close()