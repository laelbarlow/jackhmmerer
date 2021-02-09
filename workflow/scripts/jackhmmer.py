#!/usr/bin/env python3
"""Script for parsing parameters and running jackhmmer.
"""
import os
import subprocess
import pandas as pd

# Make output directory for intermediate profiles.
os.mkdir(snakemake.output.int_profile_dir)

# Define prefix for intermediate profile files.
prefix = os.path.join(snakemake.output.int_profile_dir,
        os.path.basename(snakemake.output.profile).rsplit('.', 1)[0])

# Run jackhmmer.
subprocess.call(['jackhmmer',
                 '-N', str(snakemake.params.max_iterations),
                 '--incE', str(snakemake.params.max_seq_evalue),
                 '--incdomE',
                 str(snakemake.params.max_domain_evalue['MaxEvalue'][0]),
                 '--notextw', 
                 '--chkali', prefix,
                 '-A', snakemake.output.profile,
                 '-o', snakemake.output.result_file,
                 '--domtblout', snakemake.output.tab_result_file,
                 snakemake.input.query_file,
                 snakemake.input.db_file
                ])
