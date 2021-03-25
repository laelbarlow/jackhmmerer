# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule retrieve_full_length_seqs:
    """
    Parse jackhmmer output files and retrieve sequences for listed positive
    hits from the searched database FASTA file.
    """
    input:
        result_file = \
        'results/jackhmmer/{params}/profiles/{queryfasta}_profile.stockholm',
        database_file = 'results/concat_dbs/database.faa',
        database_index_file = 'results/concat_dbs/database.faa.ssi',
        genome_id_taxon_csv = 'resources/genome_ids_by_taxon.csv'

    output:
        fasta = \
        'results/retrieve_full_length_seqs/{params}/{queryfasta}_full_length_hit_seqs.faa',
        top_hit_fasta = \
        'results/retrieve_full_length_seqs/{params}/{queryfasta}_full_length_top_hit_seqs.faa',

    conda:
        '../envs/retrieve_full_length_seqs.yaml'

    script:
        '../scripts/retrieve_full_length_seqs.py'


