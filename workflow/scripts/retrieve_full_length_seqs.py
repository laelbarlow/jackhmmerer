#!/usr/bin/env python3
"""Script for fetching full-length jackhmmer hit sequences from a very large FASTA file.
"""

import os
import subprocess
from Bio import SeqIO
from get_genome_ids_by_taxon import get_genome_ids_by_taxon


def insert_query_name_from_filename(fasta):
    """Append query name from file name to sequence IDs in a given fasta
    file.
    """
    # Extract query name from input FASTA filename.
    query_name = os.path.basename(fasta).split('_')[0]

    # Define temporary file path.
    fasta_temp = fasta + '_TEMP'

    # Iterate through sequences and make a list of modified records.
    new_seqs = []
    for seq in SeqIO.parse(fasta, 'fasta'):
        seq.id = seq.id + '__' + query_name
        new_seqs.append(seq)

    # Write seq records with new IDs to temp file.
    with open(fasta_temp, 'w') as o:
        SeqIO.write(new_seqs, o, 'fasta')

    # Overwrite original file with temp file.
    os.move(fasta_temp, fasta)



# Parse CSV file with genome IDs grouped by taxonomic groups of interest.
genome_ids_by_taxon = \
get_genome_ids_by_taxon(snakemake.input.genome_id_taxon_csv)


# Generate a list of sequence IDs for positive hits (hits included in the final
# jackhmmer profile).
ids_already_used = []
genome_ids_already_used = []
with open(snakemake.input.result_file) as infh,\
open(snakemake.output.fasta, 'a') as o1, \
open(snakemake.output.top_hit_fasta, 'a') as o2:
    for seq in SeqIO.parse(infh, 'stockholm'):
        # Extract original id from jackhmmer output ID (which indicates the
        # subsequence coordinates).
        original_id = seq.id.rsplit('/', 1)[0] 

        # Extract genome ID from sequence ID.
        genome_id = original_id.split('_')[0]

        # Don't retrieve the same full-length sequence more than once (even
        # though there may be more than one domain from the same sequence in
        # the results).
        if original_id in ids_already_used:
            pass

        else:
            # Retrieve sequences with esl-sfetch from the HMMer3 package.
            subprocess.call(['esl-sfetch',
                             snakemake.input.database_file,
                             original_id
                            ], stdout=o1)

            # Add the ID to the list.
            ids_already_used.append(original_id)

            # Determine whether a hit from this genome was already added to the
            # top-hit output file.
            if genome_id in genome_ids_already_used:
                pass

            else:
                # Retrieve sequence again and write to top-hit output file.
                subprocess.call(['esl-sfetch',
                                 snakemake.input.database_file,
                                 original_id
                                ], stdout=o2)

                # Add genome ID to the list of those already used.
                genome_ids_already_used.append(genome_id)


# Insert query names in FASTA headers.
insert_query_name_from_filename(snakemake.output.fasta)
insert_query_name_from_filename(snakemake.output.top_hit_fasta)


# Write taxon-specific full-length top-hit FASTA files.
for taxon in genome_ids_by_taxon.keys():
    # Define taxon-specific output file path.
    taxon_out_fasta = snakemake.output.top_hit_fasta.rsplit('.', 1)[0] \
            + '__' + taxon + '.faa'
    # Parse inclusive output FASTA, and write taxon-specific FASTA.
    with open(snakemake.output.top_hit_fasta) as infh,\
          open(taxon_out_fasta, 'w') as o:
        taxon_specific_seqs = []
        for seq in SeqIO.parse(infh, 'fasta'):
            # Extract genome ID from sequence ID.
            genome_id = seq.id.split('_')[0]

            if genome_id in genome_ids_by_taxon[taxon]:
                taxon_specific_seqs.append(seq)

        SeqIO.write(taxon_specific_seqs, o, 'fasta')

    # Insert query names in FASTA headers.
    insert_query_name_from_filename(taxon_out_fasta)





