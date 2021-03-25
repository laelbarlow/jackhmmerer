# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule plot_workflow:
    """
    Plot the snakemake workflow defined in the Snakefile file.
    """
    output:
        pdf = 'results/workflow_diagram.pdf',
        png = 'results/workflow_diagram.png'
    shell: 
        #'snakemake --cores 1 -p --rulegraph | dot -Tpdf > {output}'
        #snakemake --cores 1 -p --filegraph | dot -Tpdf > {output.pdf} && \
        #snakemake --cores 1 -p --filegraph | dot -Tpng > {output.png}
        """
        snakemake --cores 1 -p --rulegraph | dot -Tpdf > {output.pdf} && \
        snakemake --cores 1 -p --rulegraph | dot -Tpng > {output.png}
        """


