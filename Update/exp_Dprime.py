"""
Experiment D': Basin escape dynamics.

Start from a correctly-folding sequence, apply N random mutations,
measure the fraction of trials where the MFE structure changes.
"""

import sys
import random
import csv
import time

from rna_common import (
    metropolis_sample, structure_distance, fold_mfe,
    load_eterna100_from_url, make_log, RNA,
)


if __name__ == "__main__":
    log_path, tee = make_log("exp_Dprime")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    target = next(p for p in puzzles
                  if p["title"].strip().startswith("[CloudBeta] An Arm and a Leg"))
    structure = target["structure"]

    print("=" * 70)
    print("EXPERIMENT D': Basin escape dynamics")
    print("=" * 70)
    print(f"Puzzle: {target['title']}\n")

    # Find a correct-folding sequence
    print("Finding a correct-folding sequence...")
    start_seq = None
    for seed in range(42, 62):
        r = metropolis_sample(structure, T=1.0, lambda_dist=5.0,
                              n_steps=3000, n_samples=50, burn_in=800,
                              seed=seed)
        if r["exact_match"] > 0.5 and r["samples"]:
            start_seq = r["samples"][-1]
            print(f"Found with seed {seed}")
            break

    if start_seq is None:
        print("Failed to find correct-folding sequence. Aborting.")
        tee.logfile.close()
        sys.exit(1)

    print(f"\nStarting sequence: {start_seq}")
    mfe, _ = fold_mfe(start_seq)
    print(f"Initial distance: {structure_distance(mfe, structure)}\n")

    rows = []
    for n_mut in [1, 2, 3, 5, 8, 10]:
        escaped = 0
        trials = 100
        for trial in range(trials):
            rng = random.Random(2000 + trial)
            seq = list(start_seq)
            positions = rng.sample(range(len(seq)), n_mut)
            for pos in positions:
                bases = [b for b in "AUGC" if b != seq[pos]]
                seq[pos] = rng.choice(bases)
            mutated = ''.join(seq)
            mfe_m, _ = fold_mfe(mutated)
            if structure_distance(mfe_m, structure) > 0:
                escaped += 1
        rate = escaped / trials
        rows.append({
            "n_mutations": n_mut,
            "escape_rate": round(rate, 4),
            "trials": trials,
        })
        print(f"  {n_mut:>2} mutations -> escape rate {rate:.1%}")

    with open("results/exp_Dprime_basin_escape.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote results/exp_Dprime_basin_escape.csv")
    sys.stdout = tee.terminal
    tee.logfile.close()