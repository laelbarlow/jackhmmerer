#!/usr/bin/env python3
"""Script for quantifying jackhmmer search results, and comparing between
searches with different parameters (E-value thresholds).
"""

import os
import copy
import glob
from Bio.Seq import Seq
from Bio import SeqIO

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# Define functions for generating plots.

def generate_plot(labels,
                  taxa,
                  genome_counts_by_taxon,
                  output_path):
    """Take info to plot, and generate a clustered bar chart with matplotlib.
    """
    
    label_locations = np.arange(len(labels))  # the label locations
    #width = 0.35  # the width of the bars
    width = 0.15  # the width of the bars
    
    fig, ax = plt.subplots()
    clusters = []
    tnum = 0
    for taxon, genome_counts in zip(taxa, genome_counts_by_taxon):
        tnum += 1
        x = ax.bar(label_locations - ((width/len(taxa))*2*(len(taxa) - tnum)), genome_counts, width, label=taxon)
        clusters.append(x)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Genome count')
    ax.set_title('Genome counts by taxonomic group')
    ax.set_xticks(label_locations)
    ax.set_xticklabels(labels)
    ax.legend()
    
    # Add labels to bars.
    for cluster in clusters:
        autolabel(cluster, ax)
    
    # Adjust layout.
    fig.tight_layout()

    # Write plot to PNG file.
    plt.savefig(output_path)


def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height.
    """
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')



# Make output directory.
if not os.path.isdir(snakemake.output.plot_dir):
    os.mkdir(snakemake.output.plot_dir)


# Find info to plot.

# Parse CSV file with genome IDs grouped by taxonomic groups of interest.
genome_ids_by_taxon = {}
with open(snakemake.input.genome_id_taxon_csv) as infh:
    for i in infh:
        if not i.startswith('Taxon') and not i.startswith('\n'):
            si = i.strip().split(',')
            if si[0] not in genome_ids_by_taxon.keys():
                genome_ids_by_taxon[si[0]] = []
            elif not si[1] in genome_ids_by_taxon[si[0]]:
                genome_ids_by_taxon[si[0]] = genome_ids_by_taxon[si[0]] + [si[1]]

print(genome_ids_by_taxon)


# Group FASTA files by query.
fasta_files_by_query = {}
for f in snakemake.input.top_hit_fasta:
    # Extract info about file.
    query_name = os.path.basename(f).rsplit('_', 5)[0]

    # Initiate list when necessary.
    if query_name not in fasta_files_by_query.keys():
        fasta_files_by_query[query_name] = []

    # Add to dict.
    fasta_files_by_query[query_name] = fasta_files_by_query[query_name] + [f]


# Loop over lists of files for each query.
for query_name in fasta_files_by_query.keys():

    # Define output file path.
    output_path = os.path.join(snakemake.output.plot_dir, query_name + '.png')

    # Sort file paths by E-value threshold used.
    files = sorted(fasta_files_by_query[query_name],
            key=lambda x: \
            float(os.path.basename(os.path.dirname(x)).rsplit('_',
                5)[0].split('~')[1]))

    # Generate a list of genome ID counts for the files.
    genome_counts_by_evalue = []
    evalue_strings = []
    for f in files:
        # Add E-value string to list.
        evalue_strings.append(os.path.basename(os.path.dirname(f)).rsplit('_',
                5)[0].split('~')[1])

        # Get list of genome IDs from FASTA file contents.
        genome_ids = []
        with open(f) as infh:
            for seq in SeqIO.parse(infh, 'fasta'):
                genome_id = seq.id.split('_')[0]
                genome_ids.append(genome_id)

        # Count genome IDs per taxon (the FASTA file only contains top hit, so assuming
        # that there are no redundant genome IDs).
        counts_by_taxon = []
        for taxon in sorted(genome_ids_by_taxon.keys()):
            count = 0
            for genome_id in genome_ids:
                if genome_id in genome_ids_by_taxon[taxon]:
                   count += 1 
            counts_by_taxon.append(count)

        # Add genome count list to list of lists.
        genome_counts_by_evalue.append(counts_by_taxon)

    # Turn the genome_counts_by_taxon list of lists into a list of lists with
    # one list for each taxon instead of a list for each E-value/file.
    array = np.array(genome_counts_by_evalue)
    genome_counts_by_taxon = []
    for column in array.T:
        genome_counts_by_taxon.append(column)

    # Call function to plot info about searches with this query.
    labels = evalue_strings
    taxa = sorted(genome_ids_by_taxon.keys())
    generate_plot(labels,
                  taxa,
                  genome_counts_by_taxon,
                  output_path)

        

