#!/usr/bin/env python
import os
import sys
from itertools import tee, ifilterfalse
from collections import defaultdict

# gemini imports
import GeminiQuery
from GeminiQuery import select_formatter
from gemini_constants import *
from gemini_region import add_region_to_query
from gemini_subjects import (Subject, get_subjects, get_subjects_in_family,
                             get_family_dict)
from gemini_utils import itersubclasses
from dgidb import query_dgidb

def all_samples_predicate(args):
    """ returns a predicate that returns True if, for a variant,
    the only samples that have the variant have a given phenotype
    """
    subjects = get_subjects(args).values()
    return select_subjects_predicate(subjects, args)

def family_wise_predicate(args):
    formatter = select_formatter(args)
    families = get_family_dict(args)
    gq = GeminiQuery.GeminiQuery(args.db, out_format=formatter)
    predicates = []
    for f in families.values():
        family_names = [x.name for x in f]
        subjects = get_subjects_in_family(args, f).values()
        predicates.append(select_subjects_predicate(subjects, args,
                                                    family_names))
    def predicate(row):
        return sum([p(row) for p in predicates]) >= args.min_kindreds
    return predicate

def select_subjects_predicate(subjects, args, subset=None):
    subjects = set([s.name for s in subjects])
    predicates = []
    if "all" in args.in_subject:
        predicates.append(variant_in_all_subjects(subjects))
    if "none" in args.in_subject:
        predicates.append(variant_not_in_subjects(subjects))
    if "only" in args.in_subject:
        predicates.append(variant_only_in_subjects(subjects, subset))
    if "any" in args.in_subject:
        predicates.append(variant_in_any_subject(subjects))
    def predicate(row):
        return all([p(row) for p in predicates])
    return predicate

def variant_in_any_subject(subjects):
    def predicate(row):
        return subjects.intersection(samples_with_variant(row)) != set()
    return predicate

def variant_in_all_subjects(subjects):
    def predicate(row):
        return subjects.issubset(samples_with_variant(row))
    return predicate

def variant_only_in_subjects(subjects, subset=None):
    def predicate(row):
        if subset:
            check = set(subset).intersection(samples_with_variant(row))
        else:
            check = samples_with_variant(row)
        return check and subjects.issuperset(check)
    return predicate

def variant_not_in_subjects(subjects):
    def predicate(row):
        return subjects.intersection(samples_with_variant(row)) == set()
    return predicate

def samples_with_variant(row):
    return row['variant_samples']

def queries_variants(query):
    return "variants" in query.lower()

def get_row_predicates(args):
    """
    generate a list of predicates a row must pass in order to be
    returned from a query
    """
    predicates = []
    if args.family_wise:
        predicates.append(family_wise_predicate(args))
    elif args.sample_filter:
        predicates.append(all_samples_predicate(args))
    return predicates


def needs_genotypes(args):
    return (args.show_variant_samples or 
            args.family_wise or 
            args.sample_filter or 
            args.carrier_summary or
            args.show_families or 
            args.gt_filter)


def needs_gene(args):
    return (args.dgidb)


def add_required_columns_to_query(args):
    if args.region:
        add_region_to_query(args)


def run_query(args):
    predicates = get_row_predicates(args)
    add_required_columns_to_query(args)
    formatter = select_formatter(args)
    genotypes_needed = needs_genotypes(args)
    gene_needed = needs_gene(args)
    gq = GeminiQuery.GeminiQuery(args.db, out_format=formatter)
    gq.run(args.query, args.gt_filter, args.show_variant_samples,
           args.sample_delim, predicates, genotypes_needed,
           gene_needed, args.show_families)

    if args.use_header and gq.header:
        print gq.header

    if not args.dgidb:
        for row in gq:
            print row
    else:
        # collect a list of all the genes that need to be queried
        # from DGIdb
        genes = defaultdict()
        for row in gq:
            genes[row['gene']] = True

        # collect info from DGIdb
        dgidb_info = query_dgidb(genes)

        # rerun the query (the cursor is now consumed)
        gq = GeminiQuery.GeminiQuery(args.db, out_format=formatter)
        gq.run(args.query, args.gt_filter, args.show_variant_samples,
           args.sample_delim, predicates, genotypes_needed,
           gene_needed, args.show_families)

        # report the query results with DGIdb info added at the end.
        for row in gq:  
            print str(row) + "\t" + str(dgidb_info[row['gene']])


def query(parser, args):

    if (args.db is None):
        parser.print_help()

    if os.path.exists(args.db):
        run_query(args)

if __name__ == "__main__":
    main()
