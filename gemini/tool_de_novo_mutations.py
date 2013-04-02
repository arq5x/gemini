#!/usr/bin/env python
import sqlite3
import os
import sys
import collections
from copy import copy

import compression
import gemini_utils as util
from gemini_constants import *
import gemini_subjects as subjects


def get_de_novo_candidates(c, min_sample_depth=30):
    """
    Report candidate variants that meet appear to be de novo
    mutations in the child. We cannot distinguisj mutations that
    occured in the parental germline from those that occurred
    early in development in the child post-conception.
    """

    families = subjects.get_families(c)

    for family in families:

        query = "SELECT chrom, start, end, ref, alt, gene, \
                        impact, impact_severity, in_dbsnp, \
                        rs_ids, aaf_1kg_all, aaf_esp_all, \
                        clinvar_sig, clinvar_disease_name, \
                        clinvar_dbsource, gt_types, \
                        gt_depths, gts \
                 FROM variants \
                 WHERE impact_severity != 'LOW' \
                 AND num_het = 1"

        c.execute(query)
        all_query_cols = [str(tuple[0]) for tuple in c.description
                          if not tuple[0].startswith("gt")]

        family_genotype_mask = family.get_de_novo_filter()
        family_sample_gt_columns = family.get_subject_genotype_columns()
        family_sample_depth_columns = family.get_subject_depth_columns()
        family_sample_gt_labels = family.get_subject_genotype_labels()
        family_sample_dp_labels = family.get_subject_depth_labels()

        header = []
        header.append("family_id")
        for col in all_query_cols:
            header.append(col)
        for col in family_sample_gt_labels:
            header.append(col)
        for col in family_sample_dp_labels:
            header.append(col)
        yield header

        # report the resulting de_novo variants for this familiy
        for row in c:

            # unpack the genotype arrays so that we can interrogate
            # the genotypes present in each family member to conforming
            # to the genetic model being tested
            gt_types = compression.unpack_genotype_blob(row['gt_types'])
            gt_depths = compression.unpack_genotype_blob(row['gt_depths'])
            gts = compression.unpack_genotype_blob(row['gts'])

            # does the variant meet the a de novo model for this family?
            # if not, ignore.
            if not eval(family_genotype_mask):
                continue

            # make sure each sample's genotype had sufficient coverage.
            # otherwise, ignore
            insufficient_depth = False
            for col in family_sample_depth_columns:
                depth = int(eval(col))
                if depth < min_sample_depth:
                    insufficient_depth = True
                    break
            if insufficient_depth:
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

            # now report all of the depth columns
            for col in family_sample_depth_columns:
                result.append(str(eval(col)))

            yield result


def run(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        for result in get_de_novo_candidates(c, args.min_sample_depth):
            print '\t'.join(result)
