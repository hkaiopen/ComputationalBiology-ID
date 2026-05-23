"""
SBH Assembler: Lookahead Greedy Algorithm
==========================================
Implements the information-field greedy algorithm with single-step lookahead
to prevent dead-end traps (diffusion term in GL dynamics).

Reference: Liu & Huang (2026), Section 2.1
"""

import time
import random
from collections import Counter
from difflib import SequenceMatcher
from typing import List, Optional

def assemble_lookahead(kmers: List[str], k: int, start_kmer: str) -> Optional[str]:
    """
    Assemble DNA sequence from k-mer spectrum using lookahead greedy walk.
    
    Args:
        kmers: List of observed k-mers (may contain errors/repeats)
        k: k-mer length
        start_kmer: Starting k-mer (must be present in spectrum)
    
    Returns:
        Reconstructed sequence, or None if assembly fails
    """
    # 1. Frequency counting
    obs = Counter(kmers)
    unique = list(obs.keys())
    k2id = {km: i for i, km in enumerate(unique)}
    N = len(unique)
    
    # 2. Build De Bruijn graph adjacency
    adj = [[] for _ in range(N)]
    in_deg = [0] * N
    for i, a in enumerate(unique):
        for j, b in enumerate(unique):
            if a[1:] == b[:-1]:
                adj[i].append(j)
                in_deg[j] += 1
    
    # 3. Filter isolated nodes (error k-mers)
    rem = [obs[unique[i]] if (adj[i] or in_deg[i]) else 0 for i in range(N)]
    total = sum(rem)
    
    if start_kmer not in k2id:
        raise ValueError(f"Start k-mer '{start_kmer}' not in spectrum")
    
    # 4. Greedy walk with lookahead
    curr = k2id[start_kmer]
    path = [curr]
    rem[curr] -= 1
    
    for _ in range(total - 1):
        cand = [nxt for nxt in adj[curr] if rem[nxt] > 0]
        # One-step lookahead: keep only candidates with at least one successor
        safe = [nxt for nxt in cand if any(rem[nn] > 0 for nn in adj[nxt])]
        if not safe:
            safe = cand  # fallback if all are dead ends
        if not safe:
            break
        # Greedy selection: pick candidate with highest remaining fuel
        curr = max(safe, key=lambda x: rem[x])
        path.append(curr)
        rem[curr] -= 1
    
    # 5. Reconstruct full sequence
    seq = unique[path[0]]
    for i in range(1, len(path)):
        seq += unique[path[i]][-1]
    return seq


def run_demo():
    """Demonstrate algorithm on 5,000 bp genome with 2% error rate"""
    L = 5000
    k = 11
    error_rate = 0.02
    random.seed(42)
    
    bases = ['A', 'T', 'G', 'C']
    ref = ''.join(random.choices(bases, k=L))
    perfect = [ref[i:i+k] for i in range(L - k + 1)]
    
    # Inject errors
    noisy = []
    for km in perfect:
        noisy.append(km)
        if random.random() < error_rate:
            pos = random.randint(0, k-1)
            new_base = random.choice([b for b in bases if b != km[pos]])
            noisy.append(km[:pos] + new_base + km[pos+1:])
    for _ in range(int(len(perfect) * 0.01)):
        noisy.append(''.join(random.choices(bases, k=k)))
    random.shuffle(noisy)
    
    start = perfect[0]
    t0 = time.time()
    recon = assemble_lookahead(noisy, k, start)
    elapsed = time.time() - t0
    
    recon_set = set(recon[i:i+k] for i in range(len(recon)-k+1)) if recon else set()
    coverage = len(recon_set & set(perfect)) / len(set(perfect))
    match = SequenceMatcher(None, ref, recon).find_longest_match(
        0, len(ref), 0, len(recon)) if recon else None
    
    print(f"Reference length: {L} bp, k={k}, error rate: {error_rate}")
    print(f"Assembly time: {elapsed:.2f} s")
    print(f"Reconstructed length: {len(recon)} bp")
    print(f"k-mer coverage: {coverage:.2%}")
    if match:
        print(f"Longest contiguous match: {match.size} bp")


if __name__ == "__main__":
    run_demo()