#!/usr/bin/env python3
"""Script for converting from stockholm to FASTA format.
"""

import os
import copy
import glob
from Bio.Seq import Seq
from Bio import SeqIO

# Define a function for converting stockholm alignment files to FASTA.
def convert_stockholm_to_fasta(stockholm, fasta_aligned, fasta_unaligned):
    """Take a path to an alignment in stockholm format and write in FASTA
    format (aligned and unaligned) to given paths.
    """
    # Parse stockholm file with Biopython.
    records = list(SeqIO.parse(stockholm, "stockholm"))
    
    # Iterate over each 
    records_aligned = []
    records_unaligned = []
    for r in records:
        # Clear letter annotations.
        r.letter_annotations = {}
        # Append query title to ID.
        r.id = r.id.rsplit('/', 1)[0] + '__' + \
        os.path.basename(stockholm).split('_')[0]
        # Clear description.
        r.description = ''
        # Add to aligned SeqRecord list.
        records_aligned.append(copy.deepcopy(r))
        # Remove gap characters.
        r.seq = Seq(str(r.seq).replace('-', ''))
        # Add to unaligned SeqRecord list.
        records_unaligned.append(r)
        
    # Write sequences to output files.
    SeqIO.write(records_aligned, fasta_aligned, "fasta")
    SeqIO.write(records_unaligned, fasta_unaligned, "fasta")


# Convert final stockholm alignment file to FASTA format.
convert_stockholm_to_fasta(snakemake.input.alignment,
                           snakemake.output.aligned,
                           snakemake.output.unaligned
                           )
                           
# Compile a list of all intermediate input alignment files.
intermediate_alignments = \
glob.glob(os.path.join(snakemake.input.int_ali_dir, '*.sto'))

# Make output directories for intermediate sequence files.
os.mkdir(snakemake.output.int_aligned_dir)
os.mkdir(snakemake.output.int_unaligned_dir)

# Convert each intermediate alignment to FASTA format.
for f in intermediate_alignments:
    convert_stockholm_to_fasta(
        f,
        os.path.join(snakemake.output.int_aligned_dir,
            os.path.basename(f).rsplit('.', 1)[0] + \
            '.afaa'),
        os.path.join(snakemake.output.int_unaligned_dir,
            os.path.basename(f).rsplit('.', 1)[0] + \
            '.faa')
        )







