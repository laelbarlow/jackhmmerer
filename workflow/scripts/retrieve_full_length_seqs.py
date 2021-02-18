#!/usr/bin/env python3
"""Script for fetching full-length jackhmmer hit sequences from a very large FASTA file.
"""

import os
import subprocess
from Bio import SeqIO

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

