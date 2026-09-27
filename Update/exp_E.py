"""
Experiment E: Full Turner model verification.

For each short puzzle, generate 10 candidate sequences via MH sampling,
then verify with ViennaRNA whether their MFE structure matches the
target. Report strict (dist=0) and relaxed (dist<=3) rates.
"""

import sys
import csv
import time

from rna_common import (
    metropolis_sample, structure_distance, fold_mfe,
    parse_structure_features,
    load_eterna100_from_url, make_log, RNA,
)


if __name__ == "__main__":
    log_path, tee = make_log("exp_E")
    sys.stdout = tee
    print(f"Log: {log_path}")
    print(f"ViennaRNA: {RNA.__version__}\n")

    puzzles = load_eterna100_from_url()
    puzzles.sort(key=lambda p: len(p["structure"]))
    selected = [p for p in puzzles if len(p["structure"]) <= 60][:15]

    print("=" * 70)
    print("EXPERIMENT E: Turner model verification")
    print("=" * 70)
    print(f"Selected {len(selected)} short puzzles\n")

    rows = []
    for p in selected:
        structure = p["structure"]
        feats = parse_structure_features(structure)
        strict = 0
        relaxed = 0
        total = 0
        for seed in range(42, 52):   # 10 seeds
            r = metropolis_sample(structure, T=1.0, lambda_dist=5.0,
                                  n_steps=3000, n_samples=50, burn_in=800,
                                  seed=seed)
            if r["samples"]:
                s = r["samples"][-1]
                mfe, _ = fold_mfe(s)
                d = structure_distance(mfe, structure)
                total += 1
                if d == 0:
                    strict += 1
                if d <= 3:
                    relaxed += 1
        row = {
            "name": p["title"][:40],
            "length": feats["length"],
            "num_helices": feats["num_helices"],
            "strict_rate": round(strict / total, 4) if total else 0,
            "relaxed_rate": round(relaxed / total, 4) if total else 0,
            "n_samples": total,
        }
        rows.append(row)
        print(f"{p['title'][:30]:<30} L={feats['length']:>3} "
              f"strict={row['strict_rate']:.0%} "
              f"relaxed={row['relaxed_rate']:.0%}")

    with open("results/exp_E_turner_verification.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote results/exp_E_turner_verification.csv")
    sys.stdout = tee.terminal
    tee.logfile.close()