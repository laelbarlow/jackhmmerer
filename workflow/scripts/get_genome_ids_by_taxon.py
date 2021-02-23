#!/usr/bin/env python3
"""Module for functions shared in common between scripts.
"""

def get_genome_ids_by_taxon(csv):
    """Parse CSV file with genome IDs grouped by taxonomic groups of interest.
    """
    genome_ids_by_taxon = {}
    with open(csv) as infh:
        for i in infh:
            if not i.startswith('Taxon') and not i.startswith('\n'):
                si = i.strip().split(',')
                if si[0] not in genome_ids_by_taxon.keys():
                    genome_ids_by_taxon[si[0]] = []
                elif not si[1] in genome_ids_by_taxon[si[0]]:
                    genome_ids_by_taxon[si[0]] = genome_ids_by_taxon[si[0]] + [si[1]]
    return genome_ids_by_taxon

