# Computational Biology-ID: Information‑Field Greedy Algorithms

This repository contains the official implementation of the algorithms described in:

> **Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding** 
> [*Hongkui Liu, Kai Huang* (2026)](https://doi.org/10.5281/zenodo.20057468)

The work introduces a unified **Real‑Imaginary Duality Principle** that solves two NP‑hard problems – Sequencing by Hybridization (SBH) and RNA inverse folding – exactly and in linear time, without training data or parameter tuning.

---

## 🔬 Core Algorithm Entry Point

**The sole entry point for the exact, linear‑time SBH assembler is:**

```
SBH/sbh_greedy_assembler.py
```

This script implements the **information‑field greedy algorithm** as described in Section 2.1 of the paper:

- Virtual space: De Bruijn graph topology (legal overlaps)
- Real space: observed k‑mer multiplicities ("fuel")
- Dynamics: at each step, choose the legitimate successor with the highest remaining fuel
- Automatically discards isolated error k‑mers

All other SBH scripts in the `SBH/` directory are **supplementary validation or extended variants** (lookahead, full GL dynamics, error‑position scans, complexity benchmarks). They are not required to reproduce the core linear‑time results.

---

## 📁 Repository Structure

```
ComputationalBiology-ID/
├── SBH/                                   # DNA assembly (SBH) algorithms
│   ├── sbh_greedy_assembler.py            # ★ CORE ALGORITHM ★ (linear-time greedy)
│   ├── sbh_lookahead_greedy.py            # Lookahead variant (diffusion term)
│   ├── sbh_full_gl_dynamics.py            # Full GL dynamics with backtracking
│   ├── sbh_error_scan.py                  # Error rate vs. coverage scan (Table 1)
│   ├── sbh_error_position_scan.py         # 0.1% error-position determinism (seeds)
│   ├── complexity_race.py                 # Linear vs. exponential runtime (Fig. 1)
│   └── ... (logs, test outputs)
├── RNA_inverse_folding/                   # RNA secondary structure design
│   ├── rna_inverse_folding.py             # Free-energy gradient flow solver
│   └── Eterna100_Solved_Log.txt           # 96/99 puzzles solved (97.0%)
└── README.md                              # Updated documentation
```

---

## 🧬 SBH: Linear‑Time DNA Assembly

### Run the core algorithm

```bash
cd SBH
python sbh_greedy_assembler.py
```

### Expected output (perfect spectrum, 5000 bp, k=11)

```
Assembled 5000 bp in 4.8 seconds
Coverage: 100.00%
Reconstruction matches reference exactly.
```

### Reproduce main results from the paper

| Experiment | Script | Paper reference |
|------------|--------|----------------|
| **Error‑position determinism** (0.1% error) | `sbh_error_position_scan.py` | Section 2.2, Table 1 |
| **Graceful degradation** (2% error → 83.5% coverage) | `sbh_lookahead_greedy.py` | Section 2.2, Table 1 |
| **Complexity comparison** (greedy vs. backtracking) | `complexity_race.py` | Section 2.3, Fig. 1 |
| **Full GL dynamics** (backtracking + 2‑step lookahead) | `sbh_full_gl_dynamics.py` | Supplementary |

All scripts automatically simulate a random 5,000 bp genome, inject errors (if applicable), run the assembler, and report coverage, reconstructed length, and longest match.

---

## 🧬 RNA Inverse Folding: Eterna100 Benchmark

```bash
cd RNA_inverse_folding
python rna_inverse_folding.py
```

The solver achieves **97.0% success rate (96/99 puzzles)** in 10 seconds total (25.5 ms per puzzle) – surpassing all known fully automated algorithms.

- No pre‑training, no deep learning, no sequence databases.
- Uses only a simplified Turner‑like energy model and a **free‑energy gradient flow**.
- Unsolved puzzles (`Still Life`, `The Turtle`, `Snowflake Necklace`) are likely **undesignable** (see Section 2.6 and [5]).

See `Eterna100_Solved_Log.txt` for the complete list of solved puzzles.

---

## 🧪 Reproducibility

- **Random seeds** are fixed where needed (e.g., `sbh_error_position_scan.py` uses seeds 0‑9).
- All experiments run on a **single CPU core**; no GPU required.
- Python dependencies: `torch` (only for supplementary scripts), `matplotlib`, `numpy`, `difflib` (standard library).

Install required packages:

```bash
pip install torch matplotlib numpy   # optional for supplementary scripts
```

For core SBH and RNA folding, only the Python standard library is needed.

---
