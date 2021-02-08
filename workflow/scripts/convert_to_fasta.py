#!/usr/bin/env python3
"""Script for converting from stockholm to FASTA format.
"""

from Bio.Seq import Seq
from Bio import SeqIO

records = list(SeqIO.parse(snakemake.input.alignment, "stockholm"))
count = SeqIO.write(records, snakemake.output.aligned, "fasta")

records2 = []
for r in records:
    r.letter_annotations = {}
    r.seq = Seq(str(r.seq).replace('-', ''))
    records2.append(r)
    
SeqIO.write(records2, snakemake.output.unaligned, "fasta")

print("Converted %i records from stockholm to FASTA format." % count)
