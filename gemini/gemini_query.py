#!/usr/bin/env python

import os
import sys
from itertools import tee, ifilterfalse

# gemini imports
import GeminiQuery
from gemini_subjects import get_subjects


def variant_only_in_phenotype(args):
    """ returns a predicate that returns True if, for a variant,
    the only samples that have the variant have a given phenotype
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    subjects = get_subjects(gq.c).values()
    have_phenotype = [y.name for y in
                      filter(lambda x: x.phenotype == args.phenotype, subjects)]
    def predicate(row):
        return all(map(lambda x: x in have_phenotype, samples_with_variant(row)))
    return predicate

def variant_not_in_phenotype(args):
    """ returns a predicate that returns True if, for a variant,
    the only samples that have the variant don't have a given phenotype
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    subjects = get_subjects(gq.c).values()
    have_phenotype = [y.name for y in
                      filter(lambda x: x.phenotype == args.phenotype, subjects)]
    def predicate(row):
        return not any(map(lambda x: x in have_phenotype, samples_with_variant(row)))
    return predicate


def get_samples_with_variant(row):
   return row['variant_samples'].split(',')

def get_predicates(args):
    predicates = []
    if args.phenotype:
        predicates.append(variant_only_in_phenotype)
    if args.exclude_phenotype:
        predicates.append(variant_not_in_phenotype)

    return predicates

def needs_samples(args):
    return args.show_variant_samples or args.phenotype or args.format == "tped"

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
