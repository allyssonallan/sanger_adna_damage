#!/usr/bin/env python3
import sys
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

aln_file = sys.argv[1]
out_file = sys.argv[2]

alignment = AlignIO.read(aln_file, "fasta")

consensus_seq = ""
for i in range(alignment.get_alignment_length()):
    bases = [record.seq[i] for record in alignment]
    base_set = set(bases) - set("N-")
    if len(base_set) == 1:
        consensus_seq += base_set.pop()
    elif len(base_set) == 0:
        consensus_seq += "N"
    else:
        consensus_seq += "N"

consensus = SeqRecord(Seq(consensus_seq), id="consensus", description="")
SeqIO.write(consensus, out_file, "fasta")
