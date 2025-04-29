#!/usr/bin/env python3
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

ab1_file = sys.argv[1]
out_fasta = sys.argv[2]
out_filtered_fasta = sys.argv[3]
plot_file = sys.argv[4]

# Parameters
min_quality = 20  # you can tune it

record = SeqIO.read(ab1_file, "abi")

# Save the raw fasta
SeqIO.write(record, out_fasta, "fasta")

# Quality extraction
qualities = record.letter_annotations["phred_quality"]
sequence = str(record.seq)

# Plot quality
plt.figure(figsize=(12,4))
plt.plot(qualities, marker='.', linestyle='-')
plt.axhline(min_quality, color='red', linestyle='--')
plt.title(f'Quality scores for {ab1_file}')
plt.xlabel('Base position')
plt.ylabel('Phred Quality')
plt.savefig(plot_file)
plt.close()

# Filter low quality
filtered_seq = "".join([base if qual >= min_quality else "N" for base, qual in zip(sequence, qualities)])

# Save filtered
filtered_record = record[:]
filtered_record.seq = filtered_seq
SeqIO.write(filtered_record, out_filtered_fasta, "fasta")
