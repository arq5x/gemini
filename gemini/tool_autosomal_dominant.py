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

def get_auto_dominant_candidates(c, families):
    """
    Report candidate variants that meet an autosomal dominant
    inheritance model.
    """
    
    for family in families:
        
        query = "SELECT chrom, start, end, ref, alt, gene, \
                        impact, impact_severity, gt_types, gts \
                 FROM variants \
                 WHERE impact_severity != 'LOW'"

        c.execute(query)
        all_query_cols = [str(tuple[0]) for tuple in c.description \
                                            if not tuple[0].startswith("gt")]
                                  
        family_genotype_mask  = family.get_auto_dominant_filter()
        family_sample_columns = family.get_subject_columns()
        family_sample_labels = family.get_subject_labels()
        
        # print a header
        print "=========================="
        print "FAMILY:", family.family_id
        print "=========================="
        print '\t'.join(col for col in all_query_cols),
        print '\t'.join(col for col in family_sample_labels)
        
        # report the resulting auto_dom variants for this familiy
        
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
        get_auto_dominant_candidates(c, families)



