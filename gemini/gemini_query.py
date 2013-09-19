#!/usr/bin/env python

import os
import sys
from itertools import tee, ifilterfalse

# gemini imports
import GeminiQuery
from gemini_subjects import get_subjects



def variant_exclusive_to_affected(args):
    """ returns a predicate that returns True if, for a variant,
    the only samples that have the variant have an affected
    phenotype
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    subjects = get_subjects(gq.c).values()
    affected =  [y.name for y in filter(lambda x: x.affected, subjects)]
    def predicate(row):
        return all(map(lambda x: x in affected, samples_with_variant(row)))
    return predicate

def get_samples_with_variant(row):
   return row['variant_samples'].split(',')

def get_predicates(args):
    predicates = []
    if args.affected:
        predicates.append(variant_exclusive_to_affected)
    return predicates

def needs_samples(args):
    return args.show_variant_samples or args.affected or args.format == "tped"

def run_query(args):


    predicates = get_predicates(args)
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    gq.run(args.query, args.gt_filter, needs_samples(args),
           args.sample_delim, predicates)

    if args.use_header and gq.header:
        print gq.header

    for row in gq:
        print row

def query(parser, args):

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):
        run_query(args)


if __name__ == "__main__":
    main()
