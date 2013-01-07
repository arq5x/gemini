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

def get_auto_recessive_candidates(c):
    """
    Report candidate variants that meet an autosomal recessive
    inheritance model.
    """
    
    families = subjects.get_families(c)
    
    print families
    for family in families:
        
        query = "SELECT chrom, start, end, ref, alt, gene, \
                        impact, impact_severity, gt_types, gts \
                 FROM variants \
                 WHERE impact_severity != 'LOW'"

        c.execute(query)
        all_query_cols = [str(tuple[0]) for tuple in c.description \
                                            if not tuple[0].startswith("gt")]
                                  
        family_genotype_mask        = family.get_auto_recessive_filter()
        family_sample_gt_columns    = family.get_subject_genotype_columns()
        family_sample_gt_labels     = family.get_subject_genotype_labels()

        # skip this family if it cannot meet an autosomal_recessive model.
        if family_genotype_mask is None:
            continue

        # yield a header
        header = []
        header.append("family_id")
        for col in all_query_cols:
            header.append(col)
        for col in family_sample_gt_labels:
            header.append(col)
        yield header

        # yield the resulting auto_rec variants for this familiy
        for row in c:

            # unpack the genotype arrays so that we can interrogate
            # the genotypes present in each family member to conforming
            # to the genetic model being tested
            gt_types  = compression.unpack_genotype_blob(row['gt_types'])
            gts       = compression.unpack_genotype_blob(row['gts'])

            # skip if the variant doesn't meet a recessive model 
            # for this family
            if not eval(family_genotype_mask):
                continue

            result = []
            # first report all of the non-genotype columns
            result.append(str(family.family_id))
            for col in all_query_cols:
                if col == 'gt_types' or col == 'gts':
                    continue
                result.append(str(row[col]))

            # now report all of the genotype columns
            for col in family_sample_gt_columns:
                result.append(str(eval(col)))

            yield result


def run(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        for result in get_auto_recessive_candidates(c):
            print '\t'.join(result)



