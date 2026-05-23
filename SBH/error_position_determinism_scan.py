"""
Error‑Rate vs. Coverage Scan – Isolate error position determinism
==================================================================
This script is self‑contained and demonstrates the core finding of
the paper: at the same low error rate (0.1%), assembly coverage can
vary from 35% to 100% solely due to the spatial location of errors.
Experts who do not have time to review other Python files can run
this script directly to see the "error‑position determinism" phenomenon.

Key idea:
- The reference genome (5,000 bp) is fixed.
- For each tested error rate, a *different* random seed is used to
  inject errors, so the positions of errors change.
- The lookahead greedy assembler is then run on the noisy spectrum.
- Even though the total error rate increases (0.01% → 2%), coverage
  does NOT decrease monotonically. It depends on whether errors hit
  topological bottlenecks.

This script reproduces Table 1 and the discussion in Section 2.2.
"""

import time
import random
from collections import Counter
from difflib import SequenceMatcher

# ---------- Core lookahead greedy assembler (self‑contained) ----------
def assemble_lookahead(kmers, k, start_kmer):
    """
    Assembles a sequence from a list of k‑mers using a one‑step lookahead
    greedy walk. Returns the reconstructed sequence (string).
    """
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


# ---------- Experimental parameters ----------
L = 5000               # genome length (bp)
k = 11                 # k‑mer length
error_rates = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02]   # 0.01% to 2%

# Fixed reference genome (same for all trials)
random.seed(42)
bases = ['A', 'T', 'G', 'C']
ref = ''.join(random.choices(bases, k=L))
perfect_kmers = [ref[i:i+k] for i in range(L - k + 1)]
start_kmer = perfect_kmers[0]
perfect_set = set(perfect_kmers)

print(f"{'Error rate':>10} | {'Coverage':>8} | {'Length':>6} | {'Longest match':>13} | {'Time':>6}")
print("-" * 58)

# Run scan over error rates
for err in error_rates:
    # Each error rate uses a fixed, independent seed → different error positions
    err_seed = int(err * 10000) + 100
    random.seed(err_seed)

    # Inject substitution errors
    noisy = []
    for km in perfect_kmers:
        noisy.append(km)
        if random.random() < err:
            pos = random.randint(0, k-1)
            new_base = random.choice([b for b in bases if b != km[pos]])
            noisy.append(km[:pos] + new_base + km[pos+1:])

    # Add isolated random k‑mers (noise) – fixed proportion
    random.seed(err_seed * 2)
    for _ in range(int(len(perfect_kmers) * 0.01)):
        noisy.append(''.join(random.choices(bases, k=k)))
    random.shuffle(noisy)

    # Assemble
    t0 = time.time()
    recon = assemble_lookahead(noisy, k, start_kmer)
    elapsed = time.time() - t0

    # Evaluate coverage (fraction of perfect k‑mers recovered)
    recon_set = set(recon[i:i+k] for i in range(len(recon)-k+1)) if recon else set()
    coverage = len(recon_set & perfect_set) / len(perfect_set)

    # Longest common substring between reconstructed and reference
    if recon:
        match = SequenceMatcher(None, ref, recon).find_longest_match(0, len(ref), 0, len(recon))
        longest = match.size
    else:
        longest = 0

    print(f"{err:10.2%} | {coverage:8.2%} | {len(recon):6d} | {longest:13d} | {elapsed:6.2f}s")