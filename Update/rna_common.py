"""
Shared utilities for all RNA inverse folding experiments.

Key principles:
- ViennaRNA is imported at module load. If it's missing, the import
  fails loudly. No silent fallback to a simplified model.
- All energy functions use the modern ViennaRNA API:
  RNA.fold_compound(seq, None, RNA.OPTION_EVAL_ONLY).eval_structure(struct)
- Random state is created per seed for full reproducibility.
"""

import os
import math
import random

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import RNA   # loud failure if ViennaRNA not installed


# ============================================================
# Energy
# ============================================================

def evaluate_energy(seq, structure):
    """Turner free energy of `seq` under `structure` (kcal/mol).

    Uses the modern ViennaRNA API. Raises immediately if unavailable.
    No silent fallback.
    """
    fc = RNA.fold_compound(seq, None, RNA.OPTION_EVAL_ONLY)
    return fc.eval_structure(structure)


def fold_mfe(seq):
    """Return (mfe_structure, mfe_energy) using ViennaRNA."""
    return RNA.fold(seq)


def _extract_pairs(structure):
    pairs = set()
    stack = []
    for i, ch in enumerate(structure):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                pairs.add((stack.pop(), i))
    return pairs


def structure_distance(s1, s2):
    """Symmetric-difference size of base-pair sets (0 = identical)."""
    if s1 is None or s2 is None:
        return 0
    return len(_extract_pairs(s1).symmetric_difference(_extract_pairs(s2)))


def dual_energy(seq, target, lambda_dist):
    """E_target(seq) + lambda * d(MFE(seq), target)."""
    e = evaluate_energy(seq, target)
    mfe, _ = fold_mfe(seq)
    d = structure_distance(mfe, target) if mfe else 0
    return e + lambda_dist * d


# ============================================================
# Structure features
# ============================================================

def parse_structure_features(structure):
    n = len(structure)
    paired = sum(1 for ch in structure if ch in '()')
    num_helices = 0
    in_paired = False
    max_helix = 0
    curr_helix = 0
    for ch in structure:
        if ch in '()':
            curr_helix += 1
            if not in_paired:
                num_helices += 1
                in_paired = True
        else:
            max_helix = max(max_helix, curr_helix)
            curr_helix = 0
            in_paired = False
    max_helix = max(max_helix, curr_helix)
    max_loop = 0
    curr_loop = 0
    for ch in structure:
        if ch == '.':
            curr_loop += 1
            max_loop = max(max_loop, curr_loop)
        else:
            curr_loop = 0
    return {
        "length": n,
        "paired_ratio": paired / n if n > 0 else 0.0,
        "num_helices": num_helices,
        "max_loop": max_loop,
        "max_helix": max_helix,
    }


# ============================================================
# MH sampler
# ============================================================

def _find_partner(structure, pos):
    n = len(structure)
    if structure[pos] not in '()':
        return None
    depth = 0
    if structure[pos] == '(':
        for j in range(pos + 1, n):
            if structure[j] == '(':
                depth += 1
            elif structure[j] == ')':
                if depth == 0:
                    return j
                depth -= 1
    else:
        for j in range(pos - 1, -1, -1):
            if structure[j] == ')':
                depth += 1
            elif structure[j] == '(':
                if depth == 0:
                    return j
                depth -= 1
    return None


def propose_mutation(seq, structure, rng):
    """Single-site mutation biased toward complementary bases at paired positions."""
    n = len(seq)
    pos = rng.randrange(n)
    bases = ['A', 'U', 'G', 'C']
    seq_list = list(seq)
    old = seq_list[pos]
    if structure[pos] in '()':
        partner = _find_partner(structure, pos)
        if partner is not None:
            pb = seq_list[partner]
            if pb == 'A':
                candidates = ['U']
            elif pb == 'U':
                candidates = ['A', 'G']
            elif pb == 'G':
                candidates = ['C', 'U']
            elif pb == 'C':
                candidates = ['G']
            else:
                candidates = bases
        else:
            candidates = bases
    else:
        candidates = bases
    new = rng.choice(candidates)
    if new == old:
        return seq, False
    seq_list[pos] = new
    return ''.join(seq_list), True


def metropolis_sample(structure, T, lambda_dist,
                      n_steps=3000, n_samples=300, burn_in=800, seed=42):
    """Metropolis-Hastings sampler on the dual-energy landscape."""
    rng = random.Random(seed)
    n = len(structure)
    seq = ''.join(rng.choice('AUGC') for _ in range(n))
    curr_E = dual_energy(seq, structure, lambda_dist)

    samples = []
    sample_interval = max(1, (n_steps - burn_in) // n_samples)

    for step in range(n_steps):
        new_seq, changed = propose_mutation(seq, structure, rng)
        if changed:
            new_E = dual_energy(new_seq, structure, lambda_dist)
            dE = new_E - curr_E
            if dE <= 0 or rng.random() < math.exp(-dE / T):
                seq = new_seq
                curr_E = new_E
        if step >= burn_in and (step - burn_in) % sample_interval == 0:
            samples.append(seq)

    exact = 0
    dists = []
    for s in samples:
        mfe, _ = RNA.fold(s)
        d = structure_distance(mfe, structure)
        dists.append(d)
        if d == 0:
            exact += 1

    return {
        "exact_match": exact / len(samples) if samples else 0.0,
        "avg_distance": sum(dists) / len(dists) if dists else 0.0,
        "samples": samples,
    }


# ============================================================
# Eterna100 loader
# ============================================================

def load_eterna100_from_url():
    import requests
    import re
    url = ("https://raw.githubusercontent.com/jadeshi/SentRNA/"
           "master/data/test/eterna100.txt")
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    lines = resp.text.splitlines()

    puzzles = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if '(' in line or ')' in line:
            bracket_start = re.search(r'[\(\)]', line).start()
            title = line[:bracket_start].strip()
            rest = line[bracket_start:].strip()
            struct_end = re.search(r'[^\(\)\.]', rest)
            structure = rest[:struct_end.start()] if struct_end else rest
            if not title:
                j = i - 1
                while j >= 0 and not lines[j].strip():
                    j -= 1
                title = f"Puzzle_{len(puzzles)+1}" if j < 0 else lines[j].strip()
            puzzles.append({"title": title, "structure": structure})
        i += 1
    return [p for p in puzzles
            if '(' in p['structure'] and len(p['structure']) >= 10]


# ============================================================
# Logging
# ============================================================

class Tee:
    """Duplicate stdout writes to a log file."""
    def __init__(self, path):
        import sys
        self.terminal = sys.stdout
        self.logfile = open(path, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)
        self.logfile.flush()

    def flush(self):
        self.terminal.flush()
        self.logfile.flush()


def make_log(prefix):
    """Return (log_path, Tee instance). Caller assigns sys.stdout = tee."""
    import sys
    from datetime import datetime
    import os
    os.makedirs("results", exist_ok=True)
    path = f"results/{prefix}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    return path, Tee(path)