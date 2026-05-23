"""
SBH Assembler: Full GL Dynamics
================================
Implements the complete information-field dynamics with all four operations:
1. Information conservation (hard mask of legal overlaps)
2. Diffusion (lookahead check)
3. Gradient flow (greedy selection)
4. Anti-diffusion (local backtracking)
5. Nonlinear compression (2-step scoring)

Reference: Liu & Huang (2026), Section 4.1
"""

import time
import random
from collections import Counter
from difflib import SequenceMatcher
from typing import List, Optional, Tuple

def assemble_full_gl(kmers: List[str], k: int, start_kmer: str,
                      max_backtracks: int = 20) -> Optional[str]:
    """
    Full GL dynamics assembly with local backtracking and multi-step scoring.
    
    Args:
        kmers: Observed k-mer spectrum
        k: k-mer length
        start_kmer: Starting k-mer
        max_backtracks: Maximum number of local backtracking steps allowed
    
    Returns:
        Reconstructed sequence
    """
    obs = Counter(kmers)
    unique = list(obs.keys())
    k2id = {km: i for i, km in enumerate(unique)}
    N = len(unique)
    
    # Build De Bruijn graph
    adj = [[] for _ in range(N)]
    in_deg = [0] * N
    for i, a in enumerate(unique):
        for j, b in enumerate(unique):
            if a[1:] == b[:-1]:
                adj[i].append(j)
                in_deg[j] += 1
    
    # Filter isolated errors
    rem = [obs[unique[i]] if (adj[i] or in_deg[i]) else 0 for i in range(N)]
    total = sum(rem)
    
    if start_kmer not in k2id:
        raise ValueError(f"Start k-mer '{start_kmer}' not found")
    
    curr = k2id[start_kmer]
    path = [curr]
    rem[curr] -= 1
    backtrack_count = 0
    
    for _ in range(total - 1):
        candidates = [nxt for nxt in adj[curr] if rem[nxt] > 0]
        
        # Two-step lookahead scoring
        scored_candidates: List[Tuple[int, float]] = []
        for nxt in candidates:
            one_step = [nn for nn in adj[nxt] if rem[nn] > 0]
            if not one_step:
                continue
            
            two_step_count = 0
            for nn in one_step:
                two_step = [nnn for nnn in adj[nn] if rem[nnn] > 0]
                two_step_count += len(two_step)
            
            # Score combines: fuel + future reachability
            score = rem[nxt] + len(one_step) * 0.5 + two_step_count * 0.1
            scored_candidates.append((nxt, score))
        
        # Anti-diffusion: backtrack when stuck
        if not scored_candidates:
            if backtrack_count < max_backtracks and len(path) > 1:
                backtrack_steps = min(3, len(path) - 1)
                for _ in range(backtrack_steps):
                    last = path.pop()
                    rem[last] += 1
                curr = path[-1]
                backtrack_count += 1
                continue
            else:
                break
        
        # Gradient flow: pick highest-scoring candidate
        curr, _ = max(scored_candidates, key=lambda x: x[1])
        path.append(curr)
        rem[curr] -= 1
    
    # Reconstruct sequence
    seq = unique[path[0]]
    for i in range(1, len(path)):
        seq += unique[path[i]][-1]
    return seq


def run_demo():
    """Benchmark on 5,000 bp genome with 2% error"""
    L = 5000
    k = 11
    error_rate = 0.02
    random.seed(42)
    
    bases = ['A', 'T', 'G', 'C']
    ref = ''.join(random.choices(bases, k=L))
    perfect = [ref[i:i+k] for i in range(L - k + 1)]
    
    # Generate noisy spectrum
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
    recon = assemble_full_gl(noisy, k, start)
    elapsed = time.time() - t0
    
    recon_set = set(recon[i:i+k] for i in range(len(recon)-k+1)) if recon else set()
    coverage = len(recon_set & set(perfect)) / len(set(perfect))
    match = SequenceMatcher(None, ref, recon).find_longest_match(
        0, len(ref), 0, len(recon)) if recon else None
    
    print(f"Reference: {L} bp, k={k}, error={error_rate}")
    print(f"Time: {elapsed:.2f} s")
    print(f"Length: {len(recon)} bp, coverage: {coverage:.2%}")
    if match:
        print(f"Longest match: {match.size} bp")


if __name__ == "__main__":
    run_demo()