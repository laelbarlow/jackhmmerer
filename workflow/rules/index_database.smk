# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule index_database:
    """
    Use the HMMer3 package utility esl-sfetch to construct an index for fast retrieval
    of sequences from the concatenated database FASTA file.
    """
    input:
        db_file = 'results/concat_dbs/database.faa'

    output:
        index_file = 'results/concat_dbs/database.faa.ssi'

    conda:
        '../envs/jackhmmer.yaml'

    shell:
        """
        esl-sfetch --index {input.db_file}
        """


