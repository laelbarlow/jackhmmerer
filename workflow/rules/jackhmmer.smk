# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule jackhmmer:
    """
    Run the jackhmmer program from the HMMer3 package.
    """
    input:
        db_file = 'results/concat_dbs/database.faa',
        query_file = f'resources/FASTA_Queries/{"{queryfasta}"}.faa'

    output:
        result_file = f'results/jackhmmer/{paramspace1.wildcard_pattern}/text_outputs/{"{queryfasta}"}.txt',
        tab_result_file = f'results/jackhmmer/{paramspace1.wildcard_pattern}/tab_outputs/{"{queryfasta}"}.txt',
        profile = \
        f'results/jackhmmer/{paramspace1.wildcard_pattern}/profiles/{"{queryfasta}"}_profile.stockholm',
        int_profile_dir = \
        directory(f'results/jackhmmer/{paramspace1.wildcard_pattern}/profiles/{"{queryfasta}"}_profile_intermediates')

    params:
        max_iterations = 10,
        max_seq_evalue = 0.05,
        max_domain_evalue = paramspace1.instance

    conda:
        '../envs/jackhmmer.yaml'

    script:
        '../scripts/jackhmmer.py'


