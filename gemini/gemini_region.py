#!/usr/bin/env python
import re
import os
import sys

import GeminiQuery
from GeminiQuery import select_formatter

def _report_results(args, query, gq):
    # report the results of the region query
    gq.run(query, show_variant_samples=args.show_variant_samples)
    if args.use_header and gq.header:
        print gq.header

    for row in gq:
        print row


def get_region(args, gq):
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

    if args.columns is not None:
        query = "SELECT " + str(args.columns) + \
                    " FROM variants "
    else:
        query = "SELECT * FROM variants "

    query += "WHERE chrom = " + "'" + chrom + "'" + \
        " AND ((start BETWEEN " + start + " AND " + end + ")" +\
        " OR (end BETWEEN " + start + " AND " + end + "))"

    if args.filter:
        query += " AND " + args.filter

    query += " ORDER BY chrom, start"

    _report_results(args, query, gq)



def get_gene(args, gq):
    """
    Report all variants in a specific gene.
    """
    if args.columns is not None:
        query = "SELECT " + str(args.columns) + \
                    " FROM variants "
    else:
        query = "SELECT * FROM variants "

    query += "WHERE gene = " + "'" + args.gene + "' "

    if args.filter:
        query += " AND " + args.filter

    query += " ORDER BY chrom, start"

    _report_results(args, query, gq)

def add_region_to_query(args):
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

    where_clause = " chrom = " + "'" + chrom + "'" + \
        " AND ((start BETWEEN " + start + " AND " + end + ")" +\
        " OR (end BETWEEN " + start + " AND " + end + "))"

    args.query = _add_to_where_clause(args.query, where_clause)


def _add_to_where_clause(query, where_clause):
    where_index = query.lower().find("where")
    prefix = query[0:where_index]
    suffix = query[where_index + len("where"):]
    if where_index == -1:
        query += " WHERE " + where_clause
    else:
        query = "{0} WHERE ({1}) AND ({2})".format(prefix, suffix, where_clause)
    return query



def region(parser, args):

    if os.path.exists(args.db):

        formatter = select_formatter(args)
        gq = GeminiQuery.GeminiQuery(args.db, out_format=formatter)

        if args.region is not None and args.gene is not None:
            sys.exit('EXITING: Choose either --reg or --gene, not both.\n')
        elif args.region is not None:
            get_region(args, gq)
        elif args.gene is not None:
            get_gene(args, gq)
