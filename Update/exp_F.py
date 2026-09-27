"""
Experiment F: Fine temperature scan around T=1.0.

For the reference puzzle, scan T in [0.90, 1.10] at 0.02 steps,
30 seeds per temperature.
"""

import sys
import csv
import time
from multiprocessing import Pool

from rna_common import (
    metropolis_sample, load_eterna100_from_url, make_log, RNA,
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
    log_path, tee = make_log("exp_F")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    target = next(p for p in puzzles
                  if p["title"].strip().startswith("[CloudBeta] An Arm and a Leg"))
    structure = target["structure"]

    TEMPERATURES = [0.90 + 0.02 * i for i in range(11)]
    SEEDS = list(range(42, 72))    # 30 seeds
    N_WORKERS = 2

    print("=" * 70)
    print("EXPERIMENT F: Fine temperature scan")
    print("=" * 70)
    print(f"Puzzle: {target['title']}")
    print(f"Temperatures: {[round(t, 2) for t in TEMPERATURES]}\n")

    rows = []
    start = time.time()
    for T in TEMPERATURES:
        t0 = time.time()
        results = parallel_seeds(structure, T, SEEDS, N_WORKERS)
        exacts = [r[1] for r in results]
        dists = [r[2] for r in results]
        basin = sum(1 for d in dists if d < 1.0) / len(dists)

        row = {
            "T": round(T, 3),
            "mean_exact": round(sum(exacts) / len(exacts), 4),
            "basin_fraction": round(basin, 4),
            "mean_dist": round(sum(dists) / len(dists), 2),
        }
        rows.append(row)
        print(f"T={T:.2f}  exact={row['mean_exact']:.2%}  "
              f"basin={basin:.1%}  dist={row['mean_dist']:.2f}  "
              f"({time.time()-t0:.0f}s)")

    with open("results/exp_F_fine_temperature_scan.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["T", "mean_exact",
                                          "basin_fraction", "mean_dist"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote results/exp_F_fine_temperature_scan.csv")
    print(f"Total: {time.time()-start:.0f}s")
    sys.stdout = tee.terminal
    tee.logfile.close()