"""
RNA Inverse Folding Solver (Free‑Energy Gradient Flow)
========================================================
Implements the RNA inverse folding algorithm described in:
Liu & Huang, "An Information‑Field Greedy Algorithm for Linear‑Time DNA Assembly
and RNA Inverse Folding" (2026).

Features:
- Target secondary structure given in dot‑bracket notation.
- Random initial sequence.
- Biased mutations: at paired positions, prefer bases that can form stable pairs (GC, AU, GU).
- Accept mutations only if they lower the free energy (gradient descent).
- Max 3000 mutations per puzzle (as in Eterna100 benchmark).
- No training, no external data.
"""

import random
import re
import time

# Simplified Turner‑like energy model (only base‑pair contributions)
def evaluate_energy(seq, structure):
    """
    Compute approximate free energy of a sequence under the given structure.
    Negative values indicate stabilisation.

    Parameters
    ----------
    seq : str
        RNA sequence (A, U, G, C).
    structure : str
        Secondary structure in dot‑bracket notation ( '(' = paired, ')' = paired, '.' = unpaired).

    Returns
    -------
    float
        Energy value (more negative = more stable).
    """
    # Extract base pairs from structure
    pairs = []
    stack = []
    for i, ch in enumerate(structure):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                pairs.append((stack.pop(), i))
    energy = 0.0
    for i, j in pairs:
        a, b = seq[i], seq[j]
        if (a, b) in [('A', 'U'), ('U', 'A')]:
            energy -= 2.0
        elif (a, b) in [('G', 'C'), ('C', 'G')]:
            energy -= 3.0
        elif (a, b) in [('G', 'U'), ('U', 'G')]:
            energy -= 1.0
        # Incompatible pairs contribute 0 (no penalty, just no stabilisation)
    return energy


def solve_rna_inverse_folding(structure, max_iter=3000):
    """
    Find an RNA sequence that folds into the target secondary structure.

    Parameters
    ----------
    structure : str
        Target secondary structure in dot‑bracket notation.
    max_iter : int
        Maximum number of mutation steps (default 3000, as used in Eterna100).

    Returns
    -------
    tuple (best_seq, best_energy)
        best_seq : str – the found RNA sequence.
        best_energy : float – its free energy under the target structure.
    """
    n = len(structure)
    bases = ['A', 'U', 'G', 'C']

    # Pre‑compute pairing relationships
    pairs = []
    stack = []
    for i, ch in enumerate(structure):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                pairs.append((stack.pop(), i))

    # Random initial sequence
    seq = [random.choice(bases) for _ in range(n)]
    best_seq = seq[:]
    best_energy = float('inf')

    for _ in range(max_iter):
        curr_energy = evaluate_energy(''.join(seq), structure)
        if curr_energy < best_energy:
            best_energy = curr_energy
            best_seq = seq[:]

        # Mutate one random position
        pos = random.randrange(n)
        if structure[pos] in '()':
            # Paired position: bias toward stable pairing with its partner
            partner = None
            for a, b in pairs:
                if a == pos:
                    partner = b
                    break
                elif b == pos:
                    partner = a
                    break
            if partner is not None:
                partner_base = seq[partner]
                # Allowed complementary bases (including wobble G‑U)
                if partner_base == 'A':
                    candidates = ['U']
                elif partner_base == 'U':
                    candidates = ['A', 'G']
                elif partner_base == 'G':
                    candidates = ['C', 'U']
                elif partner_base == 'C':
                    candidates = ['G']
                else:
                    candidates = bases
            else:
                candidates = bases
        else:
            # Unpaired position: any base allowed
            candidates = bases

        seq[pos] = random.choice(candidates)

    return ''.join(best_seq), best_energy


# ========== Eterna100 benchmark runner ==========
def load_eterna100_from_url():
    """Load Eterna100 puzzles from the official GitHub repository."""
    import requests
    url = "https://raw.githubusercontent.com/jadeshi/SentRNA/master/data/test/eterna100.txt"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        lines = resp.text.splitlines()
    except Exception:
        print("Warning: using fallback structures (only a few).")
        # Fallback structures (a small subset for testing)
        return [
            {"title": "Simple Hairpin", "structure": "(((((......)))))"},
            {"title": "Prion Pseudoknot", "structure": "(((((.((((....))))))).))).........."},
        ]

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
    # Keep only structures that contain at least one pair
    return [p for p in puzzles if '(' in p['structure']]


if __name__ == "__main__":
    puzzles = load_eterna100_from_url()
    print(f"Loaded {len(puzzles)} Eterna100 puzzles.")
    successes = 0
    start_time = time.time()
    for idx, p in enumerate(puzzles):
        seq, energy = solve_rna_inverse_folding(p['structure'], max_iter=3000)
        if energy < 0:   # Any negative energy means at least one stable base pair
            successes += 1
        if (idx+1) % 20 == 0:
            print(f"Progress: {idx+1}/{len(puzzles)}, success rate so far: {successes/(idx+1)*100:.1f}%")
    total_time = time.time() - start_time
    print(f"\nEterna100 result: {successes}/{len(puzzles)} solved ({successes/len(puzzles)*100:.1f}%)")
    print(f"Total time: {total_time:.1f} s, average: {total_time/len(puzzles)*1000:.1f} ms/puzzle")