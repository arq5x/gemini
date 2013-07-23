#!/usr/bin/env python
import sqlite3
import os
import sys
import GeminiQuery
from gemini_constants import *
import gemini_subjects as subjects


def get_auto_recessive_candidates(args):
    """
    Report candidate variants that meet an autosomal recessive
    inheritance model.
    """
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    
    if args.columns is not None:
        # the user only wants to report a subset of the columns
        query = "SELECT " + str(args.columns) + " FROM variants"
    else:
        # report the kitchen sink
        query = "SELECT *" + \
                ", gts, gt_types, gt_phases, gt_depths, \
                gt_ref_depths, gt_alt_depths, gt_quals" + \
                " FROM variants"

    # add any non-genotype column limits to the where clause
    if args.filter:
        query += " WHERE " + args.filter

    # collect family info
    families = subjects.get_families(gq.c)
    family_ids = []
    family_masks = []
    family_sample_gt_labels = []
    family_sample_gt_columns = []
    for family in families:
        family_masks.append(family.get_auto_recessive_filter())
        family_sample_gt_labels.append(family.get_subject_genotype_labels())
        family_sample_gt_columns.append(family.get_subject_genotype_columns())
        family_ids.append(family.family_id)

    # run the query applying any genotype filters provided by the user.
    gq.run(query, args.gt_filter)

    # print a header
    print "family_id\tfamily_members",
    print gq.header

    # yield the resulting variants for this familiy
    for row in gq:
        
        # interrogate the genotypes present in each family member to conforming
        # to the genetic model being tested
        gt_types = row['gt_types']
        gts = row['gts']

        # test the variant for each family in the db
        for idx, fam_id in enumerate(family_ids):
            family_genotype_mask = family_masks[idx]
            family_sample_gt_label = family_sample_gt_labels[idx]
            family_sample_gt_cols = family_sample_gt_columns[idx]

            # skip if the variant doesn't meet a recessive model
            # for this family
            if not eval(family_genotype_mask):
                continue

            print str(fam_id) + "\t" + \
               ",".join([str(s) for s in family_sample_gt_label]) + "\t", \
               ",".join([str(eval(s)) for s in family_sample_gt_cols]) + "\t",
            print row

def run(parser, args):
    if os.path.exists(args.db):
        get_auto_recessive_candidates(args)
