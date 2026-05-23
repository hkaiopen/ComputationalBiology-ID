"""
Diploid Haplotype Phasing
=============================================
- Uses lookahead greedy assembler
- Proper seed selection (indegree zero) for both rounds
- Reliable evaluation: LCS length + global similarity + matching fragment preview
"""

import random
from collections import Counter
from difflib import SequenceMatcher

# ---------- Core lookahead assembler ----------
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
    """Choose a k‑mer with zero indegree (left end) if possible, else most frequent."""
    obs = Counter(kmers)
    unique = list(obs.keys())
    if not unique:
        return None
    k2id = {km: i for i, km in enumerate(unique)}
    indeg = [0] * len(unique)
    for i, a in enumerate(unique):
        for j, b in enumerate(unique):
            if a[1:] == b[:-1]:
                indeg[j] += 1
    # find first indegree zero that also has outgoing edge
    for i, deg in enumerate(indeg):
        if deg == 0:
            if any(unique[i][1:] == b[:-1] for b in unique):
                return unique[i]
    # fallback: most frequent
    return max(obs, key=obs.get)

def evaluate(seq, hapA, hapB):
    """Return detailed evaluation: LCS, ratio, and which haplotype it matches better."""
    if not seq:
        return None, None, None, None
    smA = SequenceMatcher(None, seq, hapA)
    smB = SequenceMatcher(None, seq, hapB)
    lcsA = smA.find_longest_match().size
    lcsB = smB.find_longest_match().size
    ratioA = smA.ratio()
    ratioB = smB.ratio()
    best = "A" if lcsA > lcsB else "B" if lcsB > lcsA else "tie"
    return best, max(lcsA, lcsB), ratioA, ratioB

# ---------- Haplotype generation (unchanged) ----------
def generate_haplotypes(L=1000, k=21, snp_rate=0.005, indel_rate=0.002):
    bases = ['A', 'T', 'G', 'C']
    ref = ''.join(random.choices(bases, k=L))
    hapA = list(ref)
    hapB = list(ref)
    for i in range(L):
        if random.random() < snp_rate:
            hapA[i] = random.choice([b for b in bases if b != ref[i]])
        if random.random() < snp_rate:
            hapB[i] = random.choice([b for b in bases if b != ref[i]])
    # simple indel deletion (only for demonstration)
    i = 0
    while i < len(hapA):
        if random.random() < indel_rate:
            del hapA[i]
        else:
            i += 1
    i = 0
    while i < len(hapB):
        if random.random() < indel_rate:
            del hapB[i]
        else:
            i += 1
    return ''.join(hapA), ''.join(hapB)

# ---------- Main ----------
if __name__ == "__main__":
    random.seed(42)
    L = 1000
    k = 21
    hapA, hapB = generate_haplotypes(L, k, snp_rate=0.005, indel_rate=0.002)

    kmersA = [hapA[i:i+k] for i in range(len(hapA)-k+1)]
    kmersB = [hapB[i:i+k] for i in range(len(hapB)-k+1)]
    mixed = kmersA + kmersB
    random.shuffle(mixed)

    # ---- First round: proper seed (indegree zero) ----
    seed1 = find_start_kmer(mixed, k)
    if seed1 is None:
        print("No valid seed found")
        exit(1)
    seq1, rem1 = assemble_lookahead(mixed, k, seed1)
    print(f"First haplotype assembled length: {len(seq1)} bp")
    if seq1:
        print(f"  Preview (first 100 bp): {seq1[:100]}...")

    # Evaluate first sequence
    best1, lcs1, ratioA1, ratioB1 = evaluate(seq1, hapA, hapB)
    if best1:
        print(f"  Evaluation: best matches haplotype {best1}, LCS = {lcs1} bp, "
              f"similarity to A = {ratioA1:.3f}, to B = {ratioB1:.3f}")

    # ---- Collect remaining k‑mers ----
    unique_all = list(Counter(mixed).keys())
    remaining = []
    for i, cnt in enumerate(rem1):
        if cnt > 0:
            remaining.extend([unique_all[i]] * cnt)
    print(f"Remaining k‑mers after first round: {len(remaining)}")

    # ---- Second round ----
    if remaining:
        seed2 = find_start_kmer(remaining, k)
        if seed2 is None:
            print("No seed for second round")
            exit(1)
        seq2, _ = assemble_lookahead(remaining, k, seed2)
        print(f"Second haplotype assembled length: {len(seq2)} bp")
        if seq2:
            print(f"  Preview (first 100 bp): {seq2[:100]}...")
        best2, lcs2, ratioA2, ratioB2 = evaluate(seq2, hapA, hapB)
        if best2:
            print(f"  Evaluation: best matches haplotype {best2}, LCS = {lcs2} bp, "
                  f"similarity to A = {ratioA2:.3f}, to B = {ratioB2:.3f}")
    else:
        print("No remaining k‑mers for second round")