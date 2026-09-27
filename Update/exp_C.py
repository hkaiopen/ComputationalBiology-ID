"""
Experiment C: SBH with Boltzmann sampling.

Compare greedy walk with Boltzmann walk on the De Bruijn graph.
This experiment does NOT use ViennaRNA, so it is unaffected by the
earlier API bug. Provided here for completeness.

Note: the original local run showed Boltzmann gives no advantage over
greedy. This is a meaningful negative result (SBH bottleneck is graph
topology, not sampling rule).
"""

import sys
import math
import random
import csv
import time
from collections import Counter, defaultdict

from rna_common import make_log


K = 11
SEQ_LEN = 5000
ERROR_RATES = [0.001, 0.005, 0.01]
TEMPERATURES = [0.1, 0.5, 1.0, 2.0, 5.0]
N_SEEDS = 10


def make_spectrum(sequence, k, error_rate, seed):
    rng = random.Random(seed)
    spectrum = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if rng.random() < error_rate:
            pos = rng.randrange(k)
            bases = ['A', 'T', 'G', 'C']
            bases.remove(kmer[pos])
            kmer = kmer[:pos] + rng.choice(bases) + kmer[pos + 1:]
        spectrum.append(kmer)
    return spectrum


def build_graph(spectrum):
    graph = defaultdict(list)
    for kmer in spectrum:
        graph[kmer[:-1]].append(kmer[1:])
    return graph


def greedy_walk(graph, fuel, start):
    path = [start]
    current = start
    while True:
        succs = graph.get(current, [])
        valid = [s for s in succs if fuel.get(current + s[-1], 0) > 0]
        if not valid:
            break
        nxt = max(valid, key=lambda s: fuel.get(current + s[-1], 0))
        fuel[current + nxt[-1]] -= 1
        path.append(nxt)
        current = nxt
    return path


def boltzmann_walk(graph, fuel, start, T, rng):
    path = [start]
    current = start
    while True:
        succs = graph.get(current, [])
        valid = [s for s in succs if fuel.get(current + s[-1], 0) > 0]
        if not valid:
            break
        weights = []
        for s in valid:
            f = fuel.get(current + s[-1], 0)
            weights.append(math.exp(f / T) if T > 0 else float(f > 0))
        total = sum(weights)
        if total <= 0:
            break
        r = rng.random()
        cum = 0.0
        chosen = valid[-1]
        for s, w in zip(valid, weights):
            cum += w / total
            if r <= cum:
                chosen = s
                break
        fuel[current + chosen[-1]] -= 1
        path.append(chosen)
        current = chosen
    return path


def path_to_seq(path):
    if not path:
        return ""
    return path[0] + ''.join(node[-1] for node in path[1:])


def coverage(reconstructed, reference):
    from difflib import SequenceMatcher
    sm = SequenceMatcher(None, reconstructed, reference)
    m = sm.find_longest_match(0, len(reconstructed), 0, len(reference))
    return m.size / len(reference)


if __name__ == "__main__":
    log_path, tee = make_log("exp_C")
    sys.stdout = tee
    print(f"Log: {log_path}\n")

    rng = random.Random(0)
    reference = ''.join(rng.choice('ATGC') for _ in range(SEQ_LEN))

    rows = []
    for err in ERROR_RATES:
        print(f"\nError rate = {err*100:.1f}%")

        greedy_cov = []
        for seed in range(N_SEEDS):
            spectrum = make_spectrum(reference, K, err, seed)
            fuel = Counter(spectrum)
            graph = build_graph(spectrum)
            path = greedy_walk(graph, fuel, spectrum[0][:-1])
            greedy_cov.append(coverage(path_to_seq(path), reference))
        print(f"  Greedy: mean={sum(greedy_cov)/len(greedy_cov):.1%}")
        rows.append({"error_rate": err, "method": "greedy",
                     "T": None,
                     "mean_coverage": round(sum(greedy_cov)/len(greedy_cov), 4)})

        for T in TEMPERATURES:
            covs = []
            for seed in range(N_SEEDS):
                spectrum = make_spectrum(reference, K, err, seed)
                fuel = Counter(spectrum)
                graph = build_graph(spectrum)
                rng_b = random.Random(seed * 1000 + int(T * 100))
                path = boltzmann_walk(graph, fuel, spectrum[0][:-1], T, rng_b)
                covs.append(coverage(path_to_seq(path), reference))
            mean = sum(covs) / len(covs)
            print(f"  Boltzmann T={T:.1f}: mean={mean:.1%}")
            rows.append({"error_rate": err, "method": "boltzmann",
                         "T": T, "mean_coverage": round(mean, 4)})

    with open("results/exp_C_sbh_boltzmann.csv",
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote results/exp_C_sbh_boltzmann.csv")
    sys.stdout = tee.terminal
    tee.logfile.close()