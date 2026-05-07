# ComputationalBiology-ID

**Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20057468.svg)](https://doi.org/10.5281/zenodo.20057468)

This repository contains the core code and generated figures for the paper *"Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding"* (Liu & Huang, 2026).

## What is inside

| File | Description |
|------|-------------|
| `sbh_error_scan.py` | SBH assembly coverage vs. error rate. Scans error rates from 0% to 2% on a 5000‑bp random genome (k=11). Demonstrates **error‑position determinism**: coverage depends on where errors fall, not their total number. Outputs **fig_sbh_error_scan.pdf**. |
| `complexity_race.py` | Algorithmic‑complexity comparison: the information‑field greedy algorithm (linear time) vs. classical Eulerian backtracking (exponential time). Runs on perfect k‑mer spectra of increasing length. Includes an interpretation note explaining why greedy may fail on error‑free data but dramatically outperforms backtracking under realistic noise. Outputs **complexity_race.png**. |
| `fig_sbh_error_scan.pdf` | Coverage vs. error rate (left) and coverage distribution at 0.1% error rate (right). At 0.1% error, coverage fluctuates from 35% to 100% depending solely on random seed—the signature of error‑position determinism. |
| `complexity_race.png` | Runtime comparison: greedy stays below 5 s for 5000 bp; backtracking fails beyond 1500 bp. |

## Reproducing the results

### Requirements
- Python ≥ 3.8
- NumPy, SciPy, Matplotlib

Install dependencies:
```bash
pip install numpy scipy matplotlib
```

### Run the experiments
```bash
# Error‑position determinism (generates fig_sbh_error_scan.pdf)
python sbh_error_scan.py

# Complexity comparison (generates complexity_race.png)
python complexity_race.py
```

## Key findings at a glance

1. **Error‑position determinism**: At 0.1% error rate, coverage ranges from 35% to 100% depending solely on the spatial location of errors—not their total number. This reveals the extreme sensitivity of the information field to local topological breaks and provides a new language for quantitative modeling of mutation effects.

2. **Graceful degradation**: Even at 2% error rate, coverage remains around 83.5%.

3. **Linear‑time assembly**: The greedy algorithm runs in linear time, while classical backtracking fails beyond ~1500 bp even on perfect spectra.

---

## 9. License

This project is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License** (CC BY-NC-SA 4.0).

This license allows you to:
*   **Share** — copy and redistribute the material in any medium or format.
*   **Adapt** — remix, transform, and build upon the material.

Under the following terms:
1.  **Attribution (BY)** — You must give **appropriate credit**, provide a link to the license, and **indicate if changes were made**.
2.  **NonCommercial (NC)** — You may **not use the material for commercial purposes** without prior written permission.
3.  **ShareAlike (SA)** — If you remix, transform, or build upon the material, you **must distribute your contributions under the same license**.

To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/.

---

## Citation
If you use this code or data in your research, please cite:

**Hongkui Liu, Kai Huang.** *Validation of the Real-Imaginary Duality Principle in Core Challenges of Computational Biology: From Sequencing by Hybridization to RNA Inverse Folding*. Zenodo, 2026.  
DOI: [10.5281/zenodo.20057468](https://doi.org/10.5281/zenodo.20057468)

```bibtex
@misc{liu2026cb,
  author       = {Hongkui Liu and Kai Huang},
  title        = {Validation of the Real-Imaginary Duality Principle in Core 
                  Challenges of Computational Biology: From Sequencing by 
                  Hybridization to RNA Inverse Folding},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20057468},
}
```
