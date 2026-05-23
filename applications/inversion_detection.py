"""
Inversion Detection
============================
- Uses lookahead assembler for both reference and mutant alleles
- Proper seed selection for second round (indegree zero)
"""

import random
from collections import Counter
from difflib import SequenceMatcher

# ---------- Core lookahead assembler (same as above) ----------
def assemble_lookahead(kmers, k, start_kmer):
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
        return "", None
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
    return seq, rem


def find_start_kmer(kmers, k):
    """Find a k‑mer with zero indegree (leftmost)."""
    obs = Counter(kmers)
    unique = list(obs.keys())
    k2id = {km: i for i, km in enumerate(unique)}
    indeg = [0] * len(unique)
    for i, a in enumerate(unique):
        for j, b in enumerate(unique):
            if a[1:] == b[:-1]:
                indeg[j] += 1
    for i, deg in enumerate(indeg):
        if deg == 0:
            # ensure it has at least one outgoing edge
            if any(unique[i][1:] == b[:-1] for b in unique):
                return unique[i]
    return max(obs, key=obs.get)


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp[b] for b in reversed(seq))


def generate_genome(L, gc=0.5):
    bases = ['A', 'T', 'G', 'C']
    weights = [(1-gc)/2, (1-gc)/2, gc/2, gc/2]
    return ''.join(random.choices(bases, weights=weights, k=L))


if __name__ == "__main__":
    random.seed(42)
    L = 500
    k = 15
    ref = generate_genome(L, gc=0.5)

    # Inversion: 60 bp starting at position 50
    inv_start, inv_len = 50, 60
    inv_seg = ref[inv_start:inv_start+inv_len]
    mutated = ref[:inv_start] + reverse_complement(inv_seg) + ref[inv_start+inv_len:]

    # Mixed spectrum
    kmers_ref = [ref[i:i+k] for i in range(len(ref)-k+1)]
    kmers_mut = [mutated[i:i+k] for i in range(len(mutated)-k+1)]
    mixed = kmers_ref + kmers_mut
    random.shuffle(mixed)

    # First round: reference allele (start with a k‑mer from reference)
    start_ref = kmers_ref[0]
    seq_ref, rem_ref = assemble_lookahead(mixed, k, start_ref)
    print(f"Reference assembly length: {len(seq_ref)} bp")
    sim_ref = SequenceMatcher(None, ref, seq_ref).ratio() if seq_ref else 0
    print(f"Similarity to reference: {sim_ref:.4f}")

    # Collect remaining k‑mers
    unique_all = list(Counter(mixed).keys())
    remaining = []
    for i, cnt in enumerate(rem_ref):
        if cnt > 0:
            remaining.extend([unique_all[i]] * cnt)
    print(f"Remaining k‑mers: {len(remaining)}")

    # Second round: mutant allele – choose a seed with indegree zero
    if remaining:
        seed_mut = find_start_kmer(remaining, k)
        seq_mut, _ = assemble_lookahead(remaining, k, seed_mut)
        print(f"Mutant assembly length: {len(seq_mut)} bp")
        if seq_mut:
            sim_mut = SequenceMatcher(None, mutated, seq_mut).ratio()
            print(f"Similarity to mutated genome: {sim_mut:.4f}")
            # If similarity is still low, try also assembling from the reverse complement orientation?
            # But for a simple inversion, the mutant spectrum is symmetric, so the assembler should
            # be able to walk the entire inverted region if the seed is correct.
            # Let's also compute longest common substring length as a more robust measure.
            lcs_mut = SequenceMatcher(None, mutated, seq_mut).find_longest_match().size
            print(f"Longest common substring with mutated: {lcs_mut} bp")
    else:
        print("No remaining k‑mers – inversion likely identical to reference?")