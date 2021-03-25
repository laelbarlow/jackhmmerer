# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule convert_profiles_to_fasta:
    """
    Convert the jackhmmer output profile/alignment from stockholm to FASTA
    format.
    """
    input:
        alignment = \
        'results/jackhmmer/{params}/profiles/{queryfasta}_profile.stockholm',
        int_ali_dir = \
        'results/jackhmmer/{params}/profiles/{queryfasta}_profile_intermediates'

    output:
        aligned = \
        'results/convert_profiles_to_fasta/{params}/aligned/{queryfasta}_jackhmmer_profile.afaa',
        int_aligned_dir = \
        directory('results/convert_profiles_to_fasta/{params}/aligned/{queryfasta}_jackhmmer_profile_intermediates'),
        unaligned = \
        'results/convert_profiles_to_fasta/{params}/unaligned/{queryfasta}_jackhmmer_profile_seqs.faa',
        int_unaligned_dir = \
        directory('results/convert_profiles_to_fasta/{params}/unaligned/{queryfasta}_jackhmmer_profile_seqs_intermediates')


    conda:
        '../envs/convert_profiles_to_fasta.yaml'

    script:
        '../scripts/convert_profiles_to_fasta.py'


