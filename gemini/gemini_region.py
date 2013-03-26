#!/usr/bin/env python
import sqlite3
import re
import os
import sys
import gemini_utils as util


def get_region(c, args):
    region_regex = re.compile("(\S+):(\d+)-(\d+)")

    try:
        region = region_regex.findall(args.region)[0]
    except IndexError:
        sys.exit("Malformed region (--reg) string")

    if len(region) != 3:
        sys.exit("Malformed region (--reg) string")

    chrom = region[0]
    start = region[1]
    end = region[2]

    query = "SELECT * \
             FROM variants v \
             WHERE v.chrom = " + "'" + chrom + "'" + \
        " AND ((v.start BETWEEN " + start + " AND " + end + ")" +\
        " OR (v.end BETWEEN " + start + " AND " + end + "))" + \
        "ORDER BY chrom, start"
    c.execute(query)

    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        util.get_col_names_and_indices(c.description, ignore_gt_cols=True)

    if args.use_header:
        print args.separator.join(col for col in col_names)
    for row in c:
        print args.separator.join(str(row[i]) if row[i] is not None else "."
                                  for i in non_gt_idxs)


def get_gene(c, args):
    """
    Report all variants in a specific gene.
    """
    query = "SELECT * \
    FROM variants v \
    WHERE v.gene = " + "'" + args.gene + "' " \
        "ORDER BY chrom, start"
    c.execute(query)

    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        util.get_col_names_and_indices(c.description, ignore_gt_cols=True)

    if args.use_header:
        print args.separator.join(col for col in col_names)

    for row in c:
        print args.separator.join(str(row[i]) if row[i] is not None else "."
                                  for i in non_gt_idxs)


def region(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        if args.region is not None and args.gene is not None:
            sys.exit('Choose either --reg or --gene, not both')
        elif args.region is not None:
            get_region(c, args)
        elif args.gene is not None:
            get_gene(c, args)
