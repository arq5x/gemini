#!/usr/bin/env python

import os
import sys
from itertools import tee, ifilterfalse

# gemini imports
import GeminiQuery
from gemini_subjects import get_subjects, get_family_dict
from gemini_constants import *
from gemini_region import add_region_to_query

def family_wise_in_any_subject(args):
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    families = get_family_dict(gq.c)
    predicates = []
    for f in families.values():
        subjects = [x.name for x in f]
        predicates.append(variant_in_any_subject(subjects))
    def predicate(row):
        return sum([p(row) for p in predicates]) >= args.min_families
    return predicate

def family_wise_only_in_phenotype(args):
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    families = get_family_dict(gq.c)
    predicates = []
    for f in families.values():
        subjects = subjects_with_phenotype(f, args.phenotype)
        # if all of the subjects in a family have the phenotype,
        # skip the family for consideration
        if len(subjects) == len(f):
            predicates.append(lambda x: False)
        else:
            predicates.append(variant_in_all_subjects(subjects))
    def predicate(row):
        return sum([p(row) for p in predicates]) >= args.min_families
    return predicate

def family_wise_not_in_phenotype(args):
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    families = get_family_dict(gq.c)
    predicates = []
    for f in families.values():
        subjects = subjects_with_phenotype(f, args.exclude_phenotype)
        # if all of the subjects in a family have the phenotype,
        # skip the family for consideration
        if len(subjects) == len(f):
            predicates.append(lambda x: False)
        else:
            predicates.append(variant_not_in_subjects(subjects))
    def predicate(row):
        return sum([p(row) for p in predicates]) >= args.min_families
    return predicate

def variant_in_any_subject(subjects):
    def predicate(row):
        return subjects.intersect(samples_with_variant(row)) != set()
    return predicate

def variant_in_all_subjects(subjects):
    def predicate(row):
        return subjects.issubset(samples_with_variant(row))
    return predicate

def variant_only_in_subjects(subjects):
    def predicate(row):
        return subjects.issuperset(samples_with_variant(row))
    return predicate

def variant_not_in_subjects(subjects):
    def predicate(row):
        return subjects.intersection(samples_with_variant(row)) == set()
    return predicate

def subjects_with_phenotype(subjects, phenotype):
    f = lambda x: x.phenotype == eval(phenotype)
    return set([y.name for y in filter(f, subjects)])


def variant_only_in_phenotype(args):
    """ returns a predicate that returns True if, for a variant,
    the only samples that have the variant have a given phenotype
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)

    subjects = get_subjects(gq.c).values()
    subjects = subjects_with_phenotype(subjects, args.phenotype)
    return variant_only_in_subjects(subjects)

def variant_not_in_phenotype(args):
    """ returns a predicate that returns True if, for a variant,
    the only samples that have the variant don't have a given phenotype
    """
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    subjects = get_subjects(gq.c).values()
    subjects = subjects_with_phenotype(subjects, args.exclude_phenotype)
    return variant_not_in_subjects(subjects)


def samples_with_variant(row):
   return row['variant_samples'].split(',')

def get_predicates(args):
    predicates = []
    if args.phenotype:
        if args.family_wise:
            predicates.append(family_wise_only_in_phenotype(args))
        else:
            predicates.append(variant_only_in_phenotype(args))
    if args.exclude_phenotype:
        if args.family_wise:
            predicates.append(family_wise_not_in_phenotype(args))
        else:
            predicates.append(variant_not_in_phenotype(args))
    if args.family_wise and not (args.phenotype or args.exclude_phenotype):
        predicates.append(family_wise_in_any_subject(args))

    return predicates

def needs_samples(args):
    return (args.show_variant_samples or args.phenotype or args.format == "tped"
            or args.exclude_phenotype)

def modify_query(args):
    if args.region:
        add_region_to_query(args)


def run_query(args):

    predicates = get_predicates(args)
    modify_query(args)
    gq = GeminiQuery.GeminiQuery(args.db, out_format=args.format)
    gq.run(args.query, args.gt_filter, needs_samples(args),
           args.sample_delim, predicates)

    if args.use_header and gq.header:
        print gq.header

    for row in gq:
        print row

def uppercase_arguments(args):
    if args.phenotype:
        args.phenotype = args.phenotype.upper()
    if args.exclude_phenotype:
        args.exclude_phenotype = args.exclude_phenotype.upper()
    return args

def query(parser, args):

    uppercase_arguments(args)

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):
        run_query(args)

def partition(pred, iterable):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = tee(iterable)
    return ifilterfalse(pred, t1), filter(pred, t2)

if __name__ == "__main__":
    main()
