"""
SBH Error-Position Determinism Scan
====================================
Demonstrates that assembly quality depends on error location, not total count.
At 0.1% error rate, coverage fluctuates from 35% to 100% solely by random seed.

Reference: Liu & Huang (2026), Section 2.2, Table 1
"""

import time
import random
from collections import Counter
from difflib import SequenceMatcher
from typing import List, Set, Tuple

def assemble_lookahead(kmers: List[str], k: int, start_kmer: str) -> str:
    """Lookahead greedy assembler (core algorithm)"""
    obs = Counter(kmers)
    unique = list(obs.keys())
    k2id = {km: i for i, km in enumerate(unique)}
    N = len(unique)
    
    adj = [[] for _ in range(N)]
    in_deg = [0] * N
    for i, a in enumerate(unique):
        for j, b in enumerate(unique):
            if a[1:] == b[:-1]:
                adj[i].append(j)
                in_deg[j] += 1
    
    rem = [obs[unique[i]] if (adj[i] or in_deg[i]) else 0 for i in range(N)]
    total = sum(rem)
    
    if start_kmer not in k2id:
        return ""
    
    curr = k2id[start_kmer]
    path = [curr]
    rem[curr] -= 1
    
    for _ in range(total - 1):
        cand = [nxt for nxt in adj[curr] if rem[nxt] > 0]
        safe = [nxt for nxt in cand if any(rem[nn] > 0 for nn in adj[nxt])]
        if not safe:
            safe = cand
        if not safe:
            break
        curr = max(safe, key=lambda x: rem[x])
        path.append(curr)
        rem[curr] -= 1
    
    seq = unique[path[0]]
    for i in range(1, len(path)):
        seq += unique[path[i]][-1]
    return seq


def run_error_position_scan():
    """Scan multiple seeds at fixed error rate to demonstrate determinism"""
    L = 5000
    k = 11
    error_rate = 0.001  # 0.1%
    random.seed(42)
    
    bases = ['A', 'T', 'G', 'C']
    ref = ''.join(random.choices(bases, k=L))
    perfect_kmers = [ref[i:i+k] for i in range(L - k + 1)]
    start_km = perfect_kmers[0]
    perfect_set = set(perfect_kmers)
    
    print(f"{'Seed':>6} | {'Coverage':>8} | {'Length':>6} | {'Longest match':>13} | {'Time':>6}")
    print("-" * 58)
    
    for seed in range(10):
        # Use seed-dependent noise generation
        random.seed(seed)
        noisy = []
        for km in perfect_kmers:
            noisy.append(km)
            if random.random() < error_rate:
                pos = random.randint(0, k-1)
                new_base = random.choice([b for b in bases if b != km[pos]])
                noisy.append(km[:pos] + new_base + km[pos+1:])
        
        random.seed(seed * 2)
        for _ in range(int(len(perfect_kmers) * 0.01)):
            noisy.append(''.join(random.choices(bases, k=k)))
        random.shuffle(noisy)
        
        t0 = time.time()
        recon = assemble_lookahead(noisy, k, start_km)
        elapsed = time.time() - t0
        
        recon_set = set(recon[i:i+k] for i in range(len(recon)-k+1)) if recon else set()
        coverage = len(recon_set & perfect_set) / len(perfect_set)
        match = SequenceMatcher(None, ref, recon).find_longest_match(
            0, len(ref), 0, len(recon)) if recon else None
        
        print(f"{seed:6d} | {coverage:8.2%} | {len(recon):6d} | "
              f"{match.size if match else 0:13d} | {elapsed:6.2f}s")


if __name__ == "__main__":
    run_error_position_scan()