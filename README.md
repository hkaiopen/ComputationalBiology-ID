# Computational Biology-ID: An Information-Dynamics Framework for DNA Assembly and RNA Design

> **No pretraining, no external data.** Our constraint-space sampler ...

This repository contains the official implementation of the algorithms described in:

> **Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding**
> _Hongkui Liu, Kai Huang_ (2026)

Within the constraint-space framework, two coupling rules are implemented:

- **Greedy selection** (`SBH/`, `RNA_inverse_folding/`): at each step, accept only the direction that locally improves an energy or fuel metric. This is the original framework of the paper and yields linear-time assembly and fast inverse folding.
- **Metropolis-Hastings sampling** (`Update/`): accept or reject mutations with a Boltzmann probability. This allows the sampler to escape local minima and yields substantially better results under the full Turner model, as benchmarked on Eterna100.

Both coupling rules share the same three-component structure (constraint space, real space, coupling), and the choice between them is a design decision that depends on the problem's local landscape.

This repository contains the official implementation of the algorithms described in:

> **Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding**
> [*Hongkui Liu, Kai Huang* (2026)](https://doi.org/10.5281/zenodo.20057468)

The work introduces a unified **Real-Imaginary Duality Principle** that solves two NP-hard problems — Sequencing by Hybridization (SBH) and RNA inverse folding — within a single information-dynamics framework.

---

## ⚠️ Important: Two Energy Models in This Repository

This repository contains results under **two different energy models**.
They are not interchangeable, and the choice of model must be stated
when citing any result.

| Directory | Energy Model | ViennaRNA Dependency | Purpose |
|---|---|---|---|
| `SBH/` | None (graph algorithm) | No | SBH assembly |
| `RNA_inverse_folding/` | **Simplified pairing model** | No | Original fast solver |
| `Update/` | **Full Turner nearest-neighbor** | ViennaRNA 2.7.2 | Corrected landscape analysis + Eterna100 benchmark |
| `applications/` | Mixed | Depends on script | Supplementary |

---

## 📊 Energy Model Comparison

The two energy models differ in which physical interactions they account
for. The simplified model captures only the dominant pairing contribution;
the Turner model includes the full nearest-neighbor thermodynamics
measured experimentally for RNA duplexes.

| Energy Term | Simplified Model | Turner Model | Physical Meaning |
|---|---|---|---|
| **Base-pair formation** | ✅ AU=−2, GC=−3, GU=−1 | ✅ Full parameter set | Hydrogen bonding between paired bases |
| **Stacking (nearest-neighbor)** | ❌ | ✅ ~10 kcal/mol scale | π-π interaction between adjacent base pairs |
| **Terminal mismatch** | ❌ | ✅ | Interaction at helix ends |
| **Hairpin loop** | ❌ | ✅ Size-dependent penalty | Entropic cost of closing a loop |
| **Internal loop** | ❌ | ✅ Size + asymmetry penalty | Mismatched bases within a helix |
| **Bulge loop** | ❌ | ✅ Size-dependent penalty | Unpaired bases on one strand |
| **Multibranch loop** | ❌ | ✅ Affine model | Junction of ≥3 helices |
| **Dangling ends** | ❌ | ✅ | Unpaired bases adjacent to helix ends |
| **Sequence context** | ❌ | ✅ | Energy depends on neighboring pairs |
| **Salt correction** | ❌ | ✅ | Ionic-strength dependence |
| **Empirical origin** | Arbitrary | Experimentally measured | Turner 2004+ thermodynamic tables |

**In short**: The simplified model is a crude "pair-count" energy.
The Turner model is a proper nearest-neighbor thermodynamic model.
The two differ by roughly a factor of 2–5 in predicted stability for
typical RNAs.

---

## 📈 Quantitative Impact of Energy Model Choice

The two energy models are reported on different benchmarks. The
numbers below come from the current repository and are **not
directly comparable** — they use different success criteria and
different test sets.

| Metric | Simplified Model (`RNA_inverse_folding/`) | Full Turner Model (`Update/`) |
|---|---|---|
| Eterna100 success rate | **96/99 (97.0%)** | **24/96 (25.0%)** strict MFE |
| Reference puzzle basin fraction (T=1.0) | not reported | **31%** |
| Success criterion | `energy < 0` (permissive) | `MFE == target` (strict) |
| Reported runtime | 10 s total for Eterna100 | 98 min for full benchmark |

The two success criteria measure different things:

- **Simplified model**: counts a puzzle as "solved" if the designed
  sequence has at least one stable base pair (i.e. total free energy
  is negative). This is a permissive, energy-only criterion.
- **Turner model**: requires the minimum-free-energy (MFE) structure
  predicted by ViennaRNA to match the target structure exactly. This
  is a strict, structure-based criterion.

Because the criteria differ, the 97.0% Eterna100 rate of the
simplified model and the 25.0% strict rate of the Turner model
**should not be placed on the same scale**. They describe different
quantities.

---

## 🏆 Eterna100 Benchmark (Full Turner Model)

The benchmark script `Update/benchmark_eterna100.py` evaluates the
constraint-space Boltzmann sampler on the full Eterna100 dataset
using the ViennaRNA 2.7.2 Turner nearest-neighbor model. No pretraining,
no external data, single seed (42).

### Overall results (96 valid puzzles)

| Criterion | Solved | Rate |
|---|---|---|
| **Strict MFE** (MFE structure == target) | **24/96** | **25.0%** |
| **Relaxed** (base-pair distance ≤ 3) | **48/96** | **50.0%** |
| Mean best distance | — | 17.38 |

**Runtime**: 98.4 min (single CPU, 2 workers) → 61.5 s per puzzle.

### Layered by helix count (dominant predictor)

| Category | n | Strict MFE | Relaxed (≤3 bp) | Mean distance |
|---|---|---|---|---|
| **Simple (H≤3)** | 5 | **5/5 (100%)** | **5/5 (100%)** | 0.00 |
| **Moderate (4≤H≤7)** | 18 | **9/18 (50.0%)** | **18/18 (100%)** | 0.78 |
| **Complex (H≥8)** | 73 | **10/73 (13.7%)** | 25/73 (34.2%) | 22.66 |

**Key finding**: Helix count, not sequence length, is the dominant
predictor of success. Every moderate-complexity puzzle (4–7 helices)
is solved to within 3 base pairs.

### Layered by sequence length

| Category | n | Strict MFE | Relaxed (≤3 bp) | Mean distance |
|---|---|---|---|---|
| **Short (L≤100)** | 59 | 20/59 (33.9%) | 41/59 (69.5%) | 3.83 |
| **Medium (100<L≤200)** | 13 | 2/13 (15.4%) | 4/13 (30.8%) | 13.23 |
| **Long (L>200)** | 24 | 2/24 (8.3%) | 3/24 (12.5%) | 52.92 |

### Notable exact solves

| Puzzle | Length | Helices | Runtime |
|---|---|---|---|
| **JF1** | **387 nt** | **42** | 431 s |
| Mat - Martian 2 | 213 nt | 18 | 112 s |
| Adenine | 174 nt | 20 | 63 s |
| Tilted Russian Cross | 100 nt | 8 | 18 s |

The 387-nt JF1 puzzle is solved exactly under the strict MFE criterion
in the full ViennaRNA Turner model.

### Comparison with published methods

| Method | Strict MFE | Runtime | Pretraining |
|---|---|---|---|
| **DesiRNA (2025)** | 100/100 (100%) | 24 h | None |
| **This work (strict MFE)** | 24/96 (25%) | **98 min** | **None** |
| **This work (relaxed, ≤3 bp)** | 48/96 (50%) | **98 min** | **None** |

The comparison is not like-for-like: DesiRNA uses replica-exchange Monte
Carlo with a full energy model over 24 hours, while our solver uses
single-site Metropolis-Hastings over 98 minutes. The 14.6× speed
advantage comes with a lower absolute success rate, but with **no
pretraining, no external data, and full interpretability**.

---

## 🔬 Core Algorithm Entry Point

**The sole entry point for the exact, linear-time SBH assembler is:**

    SBH/sbh_greedy_assembler.py

This script implements the **information-field greedy algorithm** as
described in Section 2.1 of the paper:

- Virtual space: De Bruijn graph topology (legal overlaps)
- Real space: observed k-mer multiplicities ("fuel")
- Dynamics: at each step, choose the legitimate successor with the highest remaining fuel
- Automatically discards isolated error k-mers

All other SBH scripts in the `SBH/` directory are supplementary
validation or extended variants (lookahead, full GL dynamics,
error-position scans, complexity benchmarks). They are not required to
reproduce the core linear-time results. **None of them depend on
ViennaRNA**, so they are unaffected by the energy-model distinction.

---

## 📁 Repository Structure

```
ComputationalBiology-ID/
├── SBH/                                    # DNA assembly (SBH) algorithms
│   ├── sbh_greedy_assembler.py             # ★ CORE (linear-time greedy)
│   ├── sbh_lookahead_greedy.py             # Lookahead variant
│   ├── sbh_full_gl_dynamics.py             # Full GL dynamics
│   ├── sbh_error_scan.py                   # Error rate vs coverage
│   ├── sbh_error_position_scan.py          # Error-position determinism
│   ├── complexity_race.py                  # Linear vs exponential runtime
│   └── ... (logs, test outputs)
│
├── RNA_inverse_folding/                    # Original simplified-model solver
│   ├── rna_inverse_folding.py              # Free-energy gradient flow
│   │                                        # (uses simplified pairing energy)
│   └── Eterna100_Solved_Log.txt            # 96/99 puzzles solved
│
├── Update/                                 # ⭐ Full Turner model results
│   ├── run_all.sh                          # ★ One-command runner
│   ├── benchmark_eterna100.py              # ★ Eterna100 full-Turner benchmark
│   ├── rna_common.py                       # Shared utilities (ViennaRNA)
│   ├── verify_grok_alignment.py            # API equivalence check
│   ├── exp_A.py                            # Temperature stability
│   ├── exp_B.py                            # Structural dependence
│   ├── exp_Bprime.py                       # Medium-length puzzles
│   ├── exp_C.py                            # SBH Boltzmann sampling
│   ├── exp_D.py                            # Sampler robustness
│   ├── exp_Dprime.py                       # Basin escape dynamics
│   ├── exp_E.py                            # Strict vs relaxed success
│   ├── exp_F.py                            # Fine temperature scan
│   └── results/                            # CSV + full logs
│       ├── benchmark_eterna100.csv         # ★ Eterna100 benchmark
│       ├── benchmark_eterna100_*.log       #   full transcript
│       ├── exp_A_temperature_stability.csv
│       ├── exp_A_T1_per_seed.csv
│       ├── exp_B_structural_dependence.csv
│       ├── exp_Bprime_medium_length.csv
│       ├── exp_C_sbh_boltzmann.csv
│       ├── exp_D_sampler_robustness.csv
│       ├── exp_Dprime_basin_escape.csv
│       ├── exp_E_turner_verification.csv
│       ├── exp_F_fine_temperature_scan.csv
│       └── *.log
│
├── applications/                           # Supplementary applications
│
└── README.md
```

---

## 🧬 SBH: Linear-Time DNA Assembly

### Run the core algorithm

    cd SBH
    python sbh_greedy_assembler.py

### Expected output (perfect spectrum, 5000 bp, k=11)

    Assembled 5000 bp in 4.8 seconds
    Coverage: 100.00%
    Reconstruction matches reference exactly.

### Reproduce main results from the paper

| Experiment | Script | Paper reference |
|---|---|---|
| **Error-position determinism** (0.1% error) | `sbh_error_position_scan.py` | Section 2.2, Table 1 |
| **Graceful degradation** (2% error → 83.5% coverage) | `sbh_lookahead_greedy.py` | Section 2.2, Table 1 |
| **Complexity comparison** | `complexity_race.py` | Section 2.3, Fig. 1 |
| **Full GL dynamics** | `sbh_full_gl_dynamics.py` | Supplementary |

**Note**: All SBH scripts are pure graph algorithms — they do not use
ViennaRNA and are unaffected by the energy-model distinction.

---

## 🧬 RNA Inverse Folding

This repository contains two complementary RNA inverse-folding pipelines
under two energy models.

### Original solver (`RNA_inverse_folding/`) — simplified model

    cd RNA_inverse_folding
    python rna_inverse_folding.py

Achieves **97.0% success rate (96/99 puzzles)** in 10 seconds total
(25.5 ms per puzzle).

- **Energy model**: Simplified pairing-only energy
  (AU = −2, GC = −3, GU = −1 kcal/mol); no stacking, loop, bulge,
  or multibranch contributions.
- No pretraining, no deep learning, no sequence databases.
- Free-energy gradient flow.
- Unsolved puzzles (`Still Life`, `The Turtle`, `Snowflake Necklace`)
  are likely **undesignable** under this model.

### Corrected solver (`Update/`) — full ViennaRNA Turner model

Run the entire suite of experiments (A–F + Eterna100 benchmark) with a
single command:

    cd Update
    bash run_all.sh

`run_all.sh` runs, in order:

1. `verify_grok_alignment.py` — API equivalence check (~10 s)
2. `benchmark_eterna100.py` — Eterna100 full-Turner benchmark (~100 min)
3. `exp_A.py` — temperature stability (~20 min)
4. `exp_B.py` — structural dependence (~10 min)
5. `exp_Bprime.py` — medium-length puzzles (~40 min)
6. `exp_C.py` — SBH Boltzmann sampling (~3 min)
7. `exp_D.py` — sampler robustness (~20 min)
8. `exp_Dprime.py` — basin escape dynamics (~15 min)
9. `exp_E.py` — strict vs relaxed success (~10 min)
10. `exp_F.py` — fine temperature scan (~15 min)

All results are saved to `Update/results/` with timestamps and full logs.

To run individual experiments:

    cd Update
    python verify_grok_alignment.py     # API equivalence check
    python benchmark_eterna100.py       # Eterna100 benchmark
    python exp_A.py                     # temperature stability
    ...

Uses **ViennaRNA 2.7.2** with the full Turner nearest-neighbor model
(see Energy Model Comparison table above). This is the appropriate
baseline for comparing against modern RNA design tools and for making
quantitative claims about landscape structure.

---

## 📊 Key Results from `Update/` (Full Turner Model)

### Eterna100 benchmark

See the section above for full results. Summary: **24/96 (25.0%) strict
MFE, 48/96 (50.0%) relaxed (≤3 bp)**, 98.4 minutes on a single CPU.

### Temperature stability (Experiment A)

Reference puzzle: `[CloudBeta] An Arm and a Leg` (L=58, 7 helices),
100 seeds per temperature:

| T (kT) | mean_exact | basin_fraction | mean_dist |
|---|---|---|---|
| 0.50 | 36.77% | 27.0% | 2.85 |
| 0.75 | 38.08% | 26.0% | 2.88 |
| **1.00** | **36.25%** | **31.0%** | **2.85** |
| 1.50 | 26.15% | 16.0% | 3.21 |
| 2.00 | 26.77% | 18.0% | 3.31 |

**Broad temperature window** (T ∈ [0.5, 1.0]), peak at T=1.0.

### Structural dependence (Experiment B)

15 short puzzles (L ≤ 60), 100 seeds each at T=1.0 — see benchmark
results above for the full-Turner picture on the entire Eterna100 set.

### Sampler robustness (Experiment D)

| Sampler | basin_fraction |
|---|---|
| **single-site MH** | **31.0%** |
| block mutation | 0.0% |
| simulated annealing | 15.0% |

**Single-site MH is optimal.**

### Basin escape (Experiment D')

| Mutations | Escape rate |
|---|---|
| 1 | 70% |
| 2 | 93% |
| 3 | 99% |
| 5+ | 100% |

**Correctly-folding sequences are fragile.**

### SBH (Experiment C) — negative result

Boltzmann sampling gives **no advantage** over greedy in SBH. The
bottleneck is graph connectivity, not sampling rule.

---

## 🧪 Reproducibility

- Random seeds are fixed where needed.
- All experiments run on a **single CPU core** (2 workers for parallel
  runs).
- Dependencies for the core results:
  - SBH: `numpy`, `difflib` (standard library)
  - Original RNA solver: Python standard library
  - `Update/`: **ViennaRNA 2.7.2**, `numpy`, `requests`, `matplotlib`

Install for `Update/`:

    pip install ViennaRNA numpy requests matplotlib

---

## 📚 Citation

If you use this code, please cite:

> Liu, H., Huang, K. (2026). *Validation of the Real-Imaginary Duality
> Principle in Core Challenges of Computational Biology: From Sequencing
> by Hybridization to RNA Inverse Folding.* (https://doi.org/10.5281/zenodo.20057468)

---

## License

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0).

You are free to share and adapt the material under the following terms:

- **Attribution** – You must give appropriate credit, provide a link to the license, and indicate if changes were made.
- **NonCommercial** – You may not use the material for commercial purposes.
- **ShareAlike** – If you remix, transform, or build upon the material, you must distribute your contributions under the same license.

For commercial use, please contact the author email hkaiopen@foxmail.com.
