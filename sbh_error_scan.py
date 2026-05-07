"""
Figure 5 – SBH assembly coverage vs. error rate (error‑position determinism)
=============================================================================
A random 5000‑bp reference genome is fragmented into k‑mers (k=11).
Random substitution errors are introduced at rates ranging from 0 to 2%.
For each rate, several independent random seeds are used to place the
errors at different positions.  Coverage of the reconstructed sequence
is recorded.  The resulting figure (fig_sbh_error_scan.pdf) shows that
assembly quality is dominated by the spatial location of errors, not by
their total number.
"""
import random
import time
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from difflib import SequenceMatcher

# ======================== Assembly engine ========================
def assemble_lookahead(kmers, k, start_kmer):
    """Information‑field greedy assembly with one‑step lookahead."""
    obs = Counter(kmers)
    unique = list(obs.keys())
    k2id = {km: i for i, km in enumerate(unique)}
    N = len(unique)

    # Build De Bruijn graph (imaginary space)
    adj = [[] for _ in range(N)]
    in_deg = [0] * N
    for i, a in enumerate(unique):
        for j, b in enumerate(unique):
            if a[1:] == b[:-1]:
                adj[i].append(j)
                in_deg[j] += 1

    # Remove isolated error nodes
    rem = [obs[unique[i]] if (adj[i] or in_deg[i]) else 0 for i in range(N)]
    total = sum(rem)
    if start_kmer not in k2id:
        return None
    curr = k2id[start_kmer]
    path = [curr]
    rem[curr] -= 1

    for _ in range(total - 1):
        cand = [nxt for nxt in adj[curr] if rem[nxt] > 0]
        # One‑step lookahead: avoid dead ends
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

# ======================== Parameters ========================
L = 5000
k = 11
bases = ['A', 'T', 'G', 'C']
random.seed(42)                     # fixed reference genome
ref = ''.join(random.choices(bases, k=L))
perfect_kmers = [ref[i:i+k] for i in range(L - k + 1)]
start_km = perfect_kmers[0]
perfect_set = set(perfect_kmers)

error_rates = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02]
n_seeds_per_rate = 5
highlight_rate = 0.001              # 0.1% – highlight for seed sensitivity
highlight_seeds = 10

results = {}                        # {rate: [coverage, ...]}

# ======================== Scan error rates ========================
for rate in error_rates:
    coverages = []
    for seed in range(n_seeds_per_rate):
        random.seed(seed)
        noisy = []
        for km in perfect_kmers:
            noisy.append(km)
            if random.random() < rate:
                pos = random.randint(0, k-1)
                new_base = random.choice([b for b in bases if b != km[pos]])
                noisy.append(km[:pos] + new_base + km[pos+1:])
        # Add a small number of completely random k‑mers
        random.seed(seed * 2)
        for _ in range(int(len(perfect_kmers) * 0.01)):
            noisy.append(''.join(random.choices(bases, k=k)))
        random.shuffle(noisy)

        recon = assemble_lookahead(noisy, k, start_km)
        if recon:
            recon_set = set(recon[i:i+k] for i in range(len(recon)-k+1))
            cov = len(recon_set & perfect_set) / len(perfect_set)
        else:
            cov = 0.0
        coverages.append(cov)
    results[rate] = coverages

# Extra seeds for 0.1% error rate
extra_covs = []
for seed in range(10, 10 + highlight_seeds):
    random.seed(seed)
    noisy = []
    for km in perfect_kmers:
        noisy.append(km)
        if random.random() < highlight_rate:
            pos = random.randint(0, k-1)
            new_base = random.choice([b for b in bases if b != km[pos]])
            noisy.append(km[:pos] + new_base + km[pos+1:])
    random.seed(seed * 2)
    for _ in range(int(len(perfect_kmers) * 0.01)):
        noisy.append(''.join(random.choices(bases, k=k)))
    random.shuffle(noisy)
    recon = assemble_lookahead(noisy, k, start_km)
    if recon:
        recon_set = set(recon[i:i+k] for i in range(len(recon)-k+1))
        cov = len(recon_set & perfect_set) / len(perfect_set)
    else:
        cov = 0.0
    extra_covs.append(cov)
results[highlight_rate].extend(extra_covs)

# ======================== Figure ========================
plt.rcParams.update({'font.size': 11, 'font.family': 'serif'})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

# Left panel: coverage vs. error rate
rates = np.array(error_rates) * 100  # convert to percent
means = [np.mean(results[r]) for r in error_rates]
stds  = [np.std(results[r]) for r in error_rates]

ax1.errorbar(rates, means, yerr=stds, fmt='o-', capsize=4, color='navy',
             markersize=7, label='Mean coverage')
highlight_idx = error_rates.index(highlight_rate)
ax1.plot([rates[highlight_idx]] * len(results[highlight_rate]),
         results[highlight_rate], 'r.', markersize=8, alpha=0.7,
         label=f'Individual seeds at {rates[highlight_idx]:.1f}%')
ax1.set_xlabel('Error rate (%)')
ax1.set_ylabel('Coverage')
ax1.set_title('Assembly coverage vs. error rate')
ax1.legend()
ax1.grid(alpha=0.3)
ax1.set_ylim(-0.05, 1.05)

# Right panel: coverage histogram at 0.1% error rate
covs_01 = results[highlight_rate]
ax2.hist(covs_01, bins=12, color='coral', edgecolor='black', alpha=0.7)
ax2.axvline(np.mean(covs_01), color='navy', linestyle='--', linewidth=2,
            label=f'Mean = {np.mean(covs_01):.2f}')
ax2.set_xlabel('Coverage')
ax2.set_ylabel('Frequency')
ax2.set_title(f'Coverage distribution at 0.1% error rate '
              f'({len(covs_01)} seeds)')
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('fig_sbh_error_scan.pdf', dpi=300, bbox_inches='tight')
plt.show()

print(f"0.1% error rate coverage range: "
      f"min={np.min(covs_01):.2f}, max={np.max(covs_01):.2f}")