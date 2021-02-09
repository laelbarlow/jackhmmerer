#!/usr/bin/env python3
"""Script for fetching full-length jackhmmer hit sequences from a very large FASTA file.
"""

import os
import subprocess
from Bio import SeqIO

# Generate a list of sequence IDs for positive hits (hits included in the final
# jackhmmer profile).
ids_already_used = []
with open(snakemake.input.result_file) as infh,\
open(snakemake.output.fasta, 'a') as o:
    for seq in SeqIO.parse(infh, 'stockholm'):
        # Extract original id from jackhmmer output ID (which indicates the
        # subsequence coordinates).
        original_id = seq.id.rsplit('/', 1)[0] 

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
                            ], stdout=o)

            # Add the ID to the list.
            ids_already_used.append(original_id)

