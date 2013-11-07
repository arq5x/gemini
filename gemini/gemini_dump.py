#!/usr/bin/env python
import sqlite3
import numpy as np
import zlib
import re
import os
import cPickle

import gemini_utils as util
from GeminiQuery import GeminiQuery


def get_variants(c, args):
    """
    Report all columns in the variant table, except for the
    genotype vectors.
    """
    query = "SELECT * FROM variants \
             ORDER BY chrom, start"
    c.execute(query)

    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        util.get_col_names_and_indices(c.description, ignore_gt_cols=True)

    if args.use_header:
        print args.separator.join(col for col in col_names)
    for row in c:
        print args.separator.join(str(row[i]) if row[i] is not None else "." \
                                              for i in non_gt_idxs )


def get_genotypes(c, args):
    """For each variant, report each sample's genotype
       on a separate line.
    """
    idx_to_sample = util.map_indices_to_samples(c)

    query = "SELECT  v.chrom, v.start, v.end, \
                     v.ref, v.alt, \
                     v.type, v.sub_type, \
                     v.aaf, v.in_dbsnp, v.gene, \
                     v.gts \
             FROM    variants v \
             ORDER BY chrom, start"
    c.execute(query)

    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        util.get_col_names_and_indices(c.description, ignore_gt_cols=True)
    col_names.append('sample')
    col_names.append('genotype')

    if args.use_header:
        print args.separator.join(col for col in col_names)
    for row in c:
        gts = np.array(cPickle.loads(zlib.decompress(row['gts'])))
        for idx, gt in enumerate(gts):
            # xrange(len(row)-1) to avoid printing v.gts
            print args.separator.join(str(row[i]) for i in xrange(len(row)-1)),
            print args.separator.join([idx_to_sample[idx], gt])


def get_samples(c, args):
    """
    Report all of the information about the samples in the DB
    """
    query = "SELECT * FROM samples"
    c.execute(query)

    (col_names, col_idxs) = util.get_col_names_and_indices(c.description)
    if args.use_header:
        print args.separator.join(col_names)
    for row in c:
        print args.separator.join(str(row[i]) if row[i] is not None else "." \
                                              for i in xrange(len(row)) )


def tfam(args):
    """
    Report the information about the samples in the DB in TFAM format:
    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    """

    query = ("select family_id, name, paternal_id, maternal_id, "
             "sex, phenotype from samples")
    gq = GeminiQuery(args.db)
    gq.run(query)
    for row in gq:
        print " ".join(map(str, [row['family_id'], row['name'], row['paternal_id'],
                        row['maternal_id'], row['sex'], row['phenotype']]))


def dump(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        if args.variants:
            get_variants(c, args)
        elif args.genotypes:
            get_genotypes(c, args)
        elif args.samples:
            get_samples(c, args)
        elif args.tfam:
            tfam(args)


