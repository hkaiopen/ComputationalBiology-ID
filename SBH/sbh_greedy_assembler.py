"""
SBH Greedy Assembler (Information‑Field Greedy Algorithm)
==========================================================
Implements the linear‑time DNA assembly algorithm described in:
Liu & Huang, "An Information‑Field Greedy Algorithm for Linear‑Time DNA Assembly
and RNA Inverse Folding" (2026).

Key properties:
- Uses a De Bruijn graph (imaginary space) with observed k‑mer frequencies as "fuel".
- At each step, selects the legitimate successor with the highest remaining fuel.
- Isolated error k‑mers (no in‑edge and no out‑edge) are automatically excluded.
- Deterministic, no training, no parameter tuning.
"""

from collections import Counter

def assemble_sbh_greedy(kmers, k, start_kmer):
    """
    Assemble a DNA sequence from a spectrum of k‑mers using greedy fuel‑maximising walk.

    Parameters
    ----------
    kmers : list of str
        Observed k‑mer spectrum (may contain errors, repeats, noise).
    k : int
        Length of each k‑mer.
    start_kmer : str
        The k‑mer used as the starting seed (must be present in the spectrum).

    Returns
    -------
    str
        Reconstructed DNA sequence.
    """
    # 1. Frequency counts as "fuel"
    obs_counts = Counter(kmers)
    unique_kmers = list(obs_counts.keys())
    kmer_to_id = {km: i for i, km in enumerate(unique_kmers)}
    N = len(unique_kmers)

    # 2. Build De Bruijn adjacency (imaginary space: legitimate transitions)
    adj = [[] for _ in range(N)]
    in_deg = [0] * N
    for i, a in enumerate(unique_kmers):
        suffix = a[1:]
        for j, b in enumerate(unique_kmers):
            if b[:-1] == suffix:
                adj[i].append(j)
                in_deg[j] += 1

    # 3. Identify and discard isolated error k‑mers (no in‑edge and no out‑edge)
    legitimate = [True] * N
    for i in range(N):
        if not adj[i] and in_deg[i] == 0:
            legitimate[i] = False

    remaining_fuel = [obs_counts[unique_kmers[i]] if legitimate[i] else 0 for i in range(N)]
    total_steps = sum(remaining_fuel)

    # 4. Start from the seed k‑mer
    if start_kmer not in kmer_to_id or not legitimate[kmer_to_id[start_kmer]]:
        raise ValueError(f"Start k‑mer '{start_kmer}' is invalid or isolated.")
    curr = kmer_to_id[start_kmer]
    path = [curr]
    remaining_fuel[curr] -= 1

    # 5. Greedy walk: always pick the legitimate successor with the highest remaining fuel
    for _ in range(total_steps - 1):
        candidates = [nxt for nxt in adj[curr] if remaining_fuel[nxt] > 0]
        if not candidates:
            break   # No further moves (should not happen for well‑behaved spectra)
        # Choose the successor with maximum remaining fuel
        curr = max(candidates, key=lambda x: remaining_fuel[x])
        path.append(curr)
        remaining_fuel[curr] -= 1

    # 6. Reconstruct the full sequence from the path of k‑mers
    seq = unique_kmers[path[0]]
    for i in range(1, len(path)):
        seq += unique_kmers[path[i]][-1]
    return seq


# ========== Example usage (as in paper) ==========
if __name__ == "__main__":
    # Example from paper: 5,000 bp random genome, k=11
    # Here we show a small test case
    true_seq = "ATGCTAGCTATG"
    k = 3
    perfect_kmers = [true_seq[i:i+k] for i in range(len(true_seq)-k+1)]
    noisy_kmers = perfect_kmers + ['AGC', 'CTA', 'CGA']   # with errors and noise
    print("True sequence:", true_seq)
    print("Observed k‑mers:", noisy_kmers)
    reconstructed = assemble_sbh_greedy(noisy_kmers, k, start_kmer='ATG')
    print("Reconstructed :", reconstructed)
    print("Success:", reconstructed == true_seq)