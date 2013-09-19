#!/usr/bin/env python

import os
import sys
from itertools import tee, ifilterfalse

# gemini imports
import GeminiQuery
from gemini_subjects import get_subjects, get_families, get_family_dict
from GeminiQuery import flatten


def multiple_families_affected(args):
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    families = get_family_dict(gq.c)
    gq.run(args.query, show_variant_samples=True)
    if args.use_header:
        print "\family_id",
        print gq.header
    for row in gq:
        samples_with_variant = row['variant_samples'].split(',')
        affected_families = filter(lambda x: all_family_affected(families[x],
                                                                 samples_with_variant),
                                   families)
        if affected_families:
            print ",".join(affected_families),
            print row


def all_family_affected(family, samples):
    affected = filter(lambda x: x.affected, family)
    unaffected = filter(lambda x: not x.affected, family)
    all_affected = all(map(lambda x: x.name in samples, affected))
    no_unaffected = not any(map(lambda x: x.name in samples, unaffected))
    return all_affected and no_unaffected

def some_family_affected_predicate(n=1):
    def predicate(family, samples):
        affected = filter(lambda x: x.affected, family)
        unaffected = filter(lambda x: not x.affected, family)
        number_affected = sum(map(lambda x: x.name in samples, affected))
        no_unaffected = not any(map(lambda x: x.name in samples, unaffected))
        return number_affected >= n and no_unaffected
    return predicate



def family_affected(args):
    """
    output the results of a query, displaying only variants which occur
    only in affected individuals and not un unaffected individuals
    as specified by the phenotype status
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    families = get_family_dict(gq.c)
    gq.run(args.query, show_variant_samples=True)
    if args.use_header:
        print "\family_id",
        print gq.header

    if args.affected_threshold:
        predicate = some_family_affected_predicate(args.affected_threshold)
    else:
        predicate = all_family_affected

    for row in gq:
        samples_with_variant = row['variant_samples'].split(',')
        affected_families = filter(lambda x: predicate(families[x],
                                                       samples_with_variant),
                                   families)
        if len(affected_families) > args.min_families:
            print ",".join(affected_families),
            print row


def affected(args):
    """
    output the results of a query, displaying only variants which occur
    only in affected individuals and not un unaffected individuals
    as specified by the phenotype status
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    subjects = get_subjects(gq.c).values()
    affected =  [y.name for y in filter(lambda x: x.affected, subjects)]

    gq.run(args.query, show_variant_samples=True)
    for row in gq:
        samples_with_variant = row['variant_samples'].split(',')
        if not all(map(lambda x: x in affected, samples_with_variant)):
            continue
        print row


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
        elif args.family_affected:
            family_affected(args)
        else:
            run_query(args)

def partition(pred, iterable):
    t1, t2 = tee(iterable)
    return ifilterfalse(pred, t1), filter(pred, t2)

if __name__ == "__main__":
    main()
