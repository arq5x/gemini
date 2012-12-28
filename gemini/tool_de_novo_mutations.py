#!/usr/bin/env python
import sqlite3
import os
import sys
import numpy as np
import cPickle
import zlib
import collections
from copy import copy

import compression
import gemini_utils as util
from gemini_constants import *
import gemini_subjects as subjects

def get_de_novo_candidates(c, families):
    """
    Report candidate variants that meet appear to be de novo
    mutations in the child. We cannot distinguisj mutations that 
    occured in the parental germline from those that occurred
    early in development in the child post-conception.
    """
    for family in families:
        
        query = "SELECT chrom, start, end, ref, alt, gene, \
                        impact, impact_severity, gt_types, gts \
                 FROM variants \
                 WHERE impact_severity != 'LOW'"

        c.execute(query)
        all_query_cols = [str(tuple[0]) for tuple in c.description \
                                            if not tuple[0].startswith("gt")]
                                  
        family_genotype_mask  = family.get_de_novo_filter()
        family_sample_columns = family.get_subject_columns()
        family_sample_labels = family.get_subject_labels()
        
        # print a header
        print "=========================="
        print "FAMILY:", family.family_id
        print "=========================="
        print '\t'.join(col for col in all_query_cols),
        print '\t'.join(col for col in family_sample_labels)
        
        # report the resulting auto_rec variants for this familiy
        print family_genotype_mask
        for row in c:
                        
            # unpack the genotype arrays so that we can interrogate
            # the genotypes present in each family member to conforming
            # to the genetic model being tested
            gt_types  = compression.unpack_genotype_blob(row['gt_types'])
            gts       = compression.unpack_genotype_blob(row['gts'])
    
            # does the variant meet the inheritance model for this family?
            if not eval(family_genotype_mask):
                continue
        
            # first report all of the non-genotype columns
            for col in all_query_cols:
                if col == 'gt_types' or col == 'gts':
                    continue
                print str(row[col]) + '\t',
            # now report all of the genotype columns
            for col in family_sample_columns:
                print str(eval(col)) + '\t',
            print


def run(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        families = subjects.get_families(c)
        get_de_novo_candidates(c, families)



