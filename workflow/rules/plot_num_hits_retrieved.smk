# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule plot_num_hits_retrieved:
    """
    Plot the number of hits retrieved from each taxonomic group of interest.
    """
    input:
        genome_id_taxon_csv = 'resources/genome_ids_by_taxon.csv',
        top_hit_fasta = \
        expand('results/retrieve_full_length_seqs/{params}/{queryfasta}__top_hits.faa',\
            queryfasta = query_names, params=paramspace1.instance_patterns)

    output:
        plot_dir = directory('results/plot_num_hits_retrieved')

    conda:
        '../envs/plot_num_hits_retrieved.yaml'

    script:
        '../scripts/plot_num_hits_retrieved.py'


