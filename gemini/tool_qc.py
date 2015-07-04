#!/usr/bin/env python
import os
from collections import Counter

import GeminiQuery
from gemini_constants import *
import gemini_subjects as subjects
import gemini_bcolz as bc

def _get_sample_sex(args):
    gq = GeminiQuery.GeminiQuery(args.db)
    query = """SELECT name, sex
               FROM   samples
            """
    gq.run(query)
    sample_sex = {}
    for row in gq:
        if row['sex'] == '1':
            sex = 'male'
        elif row['sex'] == '2':
            sex = 'female'
        else:
            sex = 'unknown'
        sample_sex[row['name']] = sex
    return sample_sex

def _get_variant_range(args, chrom):
    "Return the starting and ending variant id for a given chromosome"
    gq = GeminiQuery.GeminiQuery(args.db)
    query = """SELECT min(variant_id) as cmin, max(variant_id) as cmax
               FROM   variants
               WHERE  chrom = '%s'
            """ % chrom
    gq.run(query)
    start, end = None, None
    for row in gq:
        start, end = row['cmin'], row['cmax']
    return start, end

def check_sex(args):
    """
    Test to see if the number and type of genotype calls on the
    X and Y chromosomes make sense given the sex of the
    individual that is stated in the samples table (PED)
    """

    # what are the reported sexes of each sample?
    sample_sex = _get_sample_sex(args)
    # where do the chrX variants start and end?
    X_start, X_end = _get_variant_range(args, 'chrX')

    bcpath = bc.get_bcolz_dir(args.db)
    print '\t'.join(['sample', 'sex', 'X_homref', 'X_het', 'X_homalt', 'het_homref_ratio'])
    for sample in sample_sex:
        path = "%s/%s/%s" % (bcpath, sample, 'gt_types')
        if os.path.exists(path):
            gt_types_carray = bc.bcolz.open(path, mode="r")
        else:
            # TODO
            pass
        # retrieve the subset of genotype types for this sample
        # that are from the X chromosome
        chrX_gt_types = gt_types_carray[X_start-1:X_end]
        # tally the frequency of each genotype type
        chrX_gt_counts = Counter(chrX_gt_types)

        het_homref_ratio = float(chrX_gt_counts[1])/float(chrX_gt_counts[0])
        print '\t'.join(str(s) for s in [sample, sample_sex[sample],
            chrX_gt_counts[0], chrX_gt_counts[1],
            chrX_gt_counts[3],
            het_homref_ratio])

def run_qc(args):
    if args.mode == "sex":
        check_sex(args)

def run(parser, args):
    if os.path.exists(args.db):
        run_qc(args)
