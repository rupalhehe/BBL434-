#!/usr/bin/env python3
import sys
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import math

K = 8
WINDOW = 500
STEP = 500

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)

def kmers_in_seq(seq, k):
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            yield kmer

def sliding_windows(seq, window, step):
    for start in range(0, len(seq) - window + 1, step):
        yield start, seq[start:start+window]

def zscore(obs, mean, std):
    if std == 0:
        return 0.0
    return (obs - mean) / std

def main():
    if len(sys.argv) != 2:
        print("Usage: python script1.py genomic.fa")
        sys.exit(1)

    fasta = sys.argv[1]
    genome = read_fasta(fasta)
    genome_len = len(genome)

    print(f"Genome length: {genome_len}")

    # Background k-mer frequencies over whole genome
    bg_counts = Counter(kmers_in_seq(genome, K))
    total_bg = sum(bg_counts.values())

    bg_freq = {k: v / total_bg for k, v in bg_counts.items()}

    # Sliding window k-mer counts
    window_positions = []
    window_counts = []

    for start, win_seq in sliding_windows(genome, WINDOW, STEP):
        counts = Counter(kmers_in_seq(win_seq, K))
        window_positions.append(start + WINDOW // 2)
        window_counts.append(counts)

    # Collect all k-mers seen
    all_kmers = set(bg_counts.keys())

    # Build matrix: kmer -> list of counts per window
    kmer_window_matrix = defaultdict(list)

    for counts in window_counts:
        total = sum(counts.values())
        for kmer in all_kmers:
            kmer_window_matrix[kmer].append(counts.get(kmer, 0))

    # Compute enrichment (z-score) per k-mer per window
    kmer_zscores = {}

    for kmer, counts in kmer_window_matrix.items():
        mean = sum(counts) / len(counts)
        var = sum((c - mean) ** 2 for c in counts) / len(counts)
        std = math.sqrt(var)
        zscores = [zscore(c, mean, std) for c in counts]
        kmer_zscores[kmer] = zscores

    # Identify most overrepresented k-mers (by max z-score)
    top_kmers = sorted(
        kmer_zscores.items(),
        key=lambda x: max(x[1]),
        reverse=True
    )[:5]

    print("Top enriched k-mers:")
    for kmer, z in top_kmers:
        print(f"{kmer}  max z-score = {max(z):.2f}")

    # Plot
    plt.figure(figsize=(12, 6))
    for kmer, z in top_kmers:
        plt.plot(window_positions, z, label=kmer)

    plt.axhline(0, linestyle="--", linewidth=0.8)
    plt.xlabel("Genome position (bp)")
    plt.ylabel("Z-score enrichment")
    plt.title(f"ORI signal checker – k={K}, window={WINDOW}, step={STEP}")
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

