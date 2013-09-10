#!/usr/bin/env python

import os
import sys

# gemini imports
import GeminiQuery
from gemini_subjects import get_subjects


def affected(args):
    """
    output the results of a query, displaying only variants which occur
    only in affected individuals and not un unaffected individuals
    as specified by the phenotype status
    """
    gq = GeminiQuery.GeminiQuery(args.db)
    subjects = get_subjects(gq.c).values()
    affected = [s.name for s in filter(lambda x: x.affected, subjects)]

    gq.run(args.query, show_variant_samples=True)
    processed = 0
    for row in gq:
        processed = 1
        if processed % 10000 == 0:
            sys.stderr.write("%d variants processed.\n" % (processed))
        if not all(map(lambda x: x in affected, row.gts)):
            continue
        print row


def tped(args):
    """
    output the results of a query in TPED format
    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    """

    VALID_CHROMOSOMES = map(str, range(1, 23)) + ["X", "Y", "XY", "MT"]
    NEED_COLUMNS = ["chrom", "rs_ids", "start", "gts"]
    gq = GeminiQuery.GeminiQuery(args.db)

    gq.run("select name from samples")
    sample_names = [row['name'] for row in gq]

    gq = GeminiQuery.GeminiQuery(args.db)
    query = gq.ensure_columns(args.query, NEED_COLUMNS)

    if args.use_header:
        print " ".join(["chrom", "name", "distance", "start",
                        " ".join(sample_names)])

    gq.run(query, show_variant_samples=True)
    for row in gq:
        chrom = row['chrom'].split("chr")[1]
        chrom = chrom if chrom in VALID_CHROMOSOMES else "0"
        name = row['rs_ids'] if row['rs_ids'] else "."
        start = str(row['start'])
        genotypes = " ".join(row['gts'])
        print " ".join([chrom, name, "0", start, genotypes])


def run_query(args):

    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    gq.run(args.query, args.gt_filter,
           args.show_variant_samples,
           args.sample_delim)

    if args.use_header and gq.header:
        print gq.header

    for row in gq:
        print row

def query(parser, args):

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):
        if args.affected:
            affected(args)
        else:
            run_query(args)

if __name__ == "__main__":
    main()
