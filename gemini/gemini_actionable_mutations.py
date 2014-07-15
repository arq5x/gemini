#!/usr/bin/env python
import sqlite3
from collections import defaultdict
from gemini_constants import *
import gemini_subjects
from dgidb import query_dgidb
import GeminiQuery

def get_actionable_mutations(parser, args):

    t_n_pairs = gemini_subjects.get_families(args.db)

    query = "SELECT variants.chrom, start, end, ref, alt, \
                    variants.gene, impact, is_somatic, \
                    gene_summary.in_cosmic_census \
             FROM variants, gene_summary \
             WHERE variants.is_somatic = 1 \
             AND variants.is_exonic = 1 \
             AND variants.chrom = gene_summary.chrom \
             AND variants.gene = gene_summary.gene \
             AND gene_summary.in_cosmic_census = 1"


    # collect the relevant genes and query DGIDB
    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(query)

    genes = defaultdict()
    for row in gq:
      genes[row['gene']] = True
    # collect info from DGIdb
    dgidb_info = query_dgidb(genes)


    # now rerun the query and report actionable mutations per DGIDB and COSMIC census.
    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(query)
    print'\t'.join(['tum_name', 'chrom', 'start', 'end', 'ref', 'alt', \
                    'gene', 'impact', 'is_somatic', 'in_cosmic_census', 'dgidb_info'])
    for row in gq:

        for pair in t_n_pairs:
            samples = pair.subjects
            if len(samples) != 2:
                continue

            tumor = pair.subjects[0]
            normal = pair.subjects[1]
            # swap if we guessed the tumor incorrectly
            if tumor.affected is False:
                tumor, normal = normal, tumor

            print'\t'.join(str(s) for s in [tumor.name, row['chrom'], \
                                            row['start'], row['end'], \
                                            row['ref'], row['alt'], \
                                            row['gene'], row['impact'], \
                                            row['is_somatic'], \
                                            row['in_cosmic_census'], \
                                            str(dgidb_info[row['gene']])])
