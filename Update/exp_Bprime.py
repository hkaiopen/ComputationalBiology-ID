"""
Experiment B': Basin structure of medium-length puzzles (L ~ 100 nt).

Select 3 puzzles closest to L=100 and run 50 seeds at T=1.0.
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
                          n_steps=5000, n_samples=400, burn_in=1500,
                          seed=seed)
    return seed, r["exact_match"], r["avg_distance"]


def parallel_seeds(structure, T, seeds, n_workers):
    with Pool(processes=n_workers) as pool:
        return pool.map(_worker, [(structure, T, s) for s in seeds])


if __name__ == "__main__":
    log_path, tee = make_log("exp_Bprime")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    puzzles.sort(key=lambda p: abs(len(p["structure"]) - 100))
    selected = puzzles[:3]

    print("=" * 70)
    print("EXPERIMENT B': Basin structure at L ~ 100 nt")
    print("=" * 70)

    SEEDS = list(range(42, 92))   # 50 seeds
    N_WORKERS = 2
    rows = []
    start = time.time()

    for p in selected:
        feats = parse_structure_features(p["structure"])
        t0 = time.time()
        print(f"\nPuzzle: {p['title'][:50]}")
        print(f"  L={feats['length']}, helices={feats['num_helices']}")

        results = parallel_seeds(p["structure"], 1.0, SEEDS, N_WORKERS)
        exacts = [r[1] for r in results]
        dists = [r[2] for r in results]
        basin = sum(1 for d in dists if d < 1.0) / len(dists)

        row = {
            "name": p["title"][:40],
            "length": feats["length"],
            "num_helices": feats["num_helices"],
            "max_loop": feats["max_loop"],
            "mean_exact": round(sum(exacts) / len(exacts), 4),
            "basin_fraction": round(basin, 4),
            "mean_dist": round(sum(dists) / len(dists), 2),
        }
        rows.append(row)
        print(f"  exact={row['mean_exact']:.2%}  "
              f"basin={basin:.1%}  dist={row['mean_dist']:.2f}  "
              f"({time.time()-t0:.0f}s)")

    with open("results/exp_Bprime_medium_length.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nTotal: {time.time()-start:.0f}s")
    sys.stdout = tee.terminal
    tee.logfile.close()