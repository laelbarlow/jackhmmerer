#!/usr/bin/env Rscript
# Original source: https://github.com/laelbarlow/jackhmmer
# License notice: MIT License, copyright (c) 2021 Lael D. Barlow

# Script for generating a group bar graph.

# Import libraries.
library(ggplot2)
library(R.devices)

# Use a null device (prevents writing of unnecessary output image files).
nulldev()

# Take command-line arguments.
args = commandArgs(trailingOnly=TRUE)

# Import data from CSV file.
data <- read.csv(file = args[1])

# Convert Evalue column to factor instead of numerical.
data$Evalue <- as.factor(data$Evalue)

# Change order of taxa in E-value groups (DATASET-SPECIFIC).
data$Taxon <- factor(data$Taxon, levels = c('Eukaryotes', 'Archaea', 'Bacteria'))
 
# Grouped Bar Graph.
ggplot(data, aes(fill=Taxon, y=Count, x=Evalue)) + 
    geom_bar(position="dodge", stat="identity", colour="black", width=0.65) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Taxonomic distribution of sequences retrieved using jackhmmer\nwith different E-value thresholds") +
    labs(y="Number of genomes yielding at least one\nretrieved sequence (count)", x = "E-value threshold") + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_brewer(palette="Blues", direction=-1)

# Write figure to given PDF file path.
ggsave(args[2])

# Turn off null device.
dev.off()

