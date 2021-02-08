#!/usr/bin/env python3
"""Script for parsing parameters and running jackhmmer.
"""
import subprocess
import pandas as pd

subprocess.call(['jackhmmer',
                 '-N', str(snakemake.params.max_iterations),
                 '--incE', str(snakemake.params.max_seq_evalue),
                 '--incdomE',
                 str(snakemake.params.max_domain_evalue['MaxEvalue'][0]),
                 '--notextw', 
                 '-A', snakemake.output.profile,
                 '-o', snakemake.output.result_file,
                 '--domtblout', snakemake.output.tab_result_file,
                 snakemake.input.query_file,
                 snakemake.input.db_file
                ])
