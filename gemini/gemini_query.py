#!/usr/bin/env python

import os
import sys

# gemini imports
import GeminiQuery


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
    sample_ids = [gq.sample_to_idx[x] for x in sample_names]

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
    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(args.query, args.gt_filter,
           args.show_variant_samples,
           args.sample_delim)

    if args.use_header and not args.use_json:
        print gq.header

    for row in gq:
        if args.use_json:
            print row.get_json()
        else:
            print row

def query(parser, args):

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):
        if args.tped:
            tped(args)
        else:
            run_query(args)

if __name__ == "__main__":
    main()
