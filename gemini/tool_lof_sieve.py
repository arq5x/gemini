#!/usr/bin/env python

import os
import sys
import sqlite3
import numpy as np
import cPickle
import zlib
from collections import defaultdict
from gemini.config import read_gemini_config
import gemini_utils as util

config = read_gemini_config()
path_dirname = config["annotation_dir"]


def get_ind_lof(c, args):

    idx_to_sample = util.map_indicies_to_samples(c)

    query = "SELECT DISTINCT v.chrom, v.start, v.end, v.ref, v.alt, \
                             v.most_severe_impact, v.aa_change, v.aa_length, \
                             v.gt_types, v.gts, i.gene, \
                             i.transcript \
             FROM variants v, variant_impacts i \
             WHERE v.variant_id = i.variant_id \
             AND i.is_lof='1'"

    c.execute(query)

    # # header
    # print '\t'.join(['chrom', 'start', 'end', 'ref', 'alt', \
    #                  'highest_impact', 'sample', 'genotype', \
    #                  'gene', 'transcript', 'pathway'])

    def _report_variant_pathways(c, args, idx_to_sample):

        (agn_paths, hgnc_paths, ensembl_paths) = get_pathways(args)

        for r in c:
            gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
            gts      = np.array(cPickle.loads(zlib.decompress(r['gts'])))        
            gene     = str(r['gene'])
            trans    = str(r['transcript'])

            for idx, type in enumerate(gt_types):
                if type > 0 and len(pathways) > 0:
                    print "\t".join([r['chrom'], str(r['start']), \
                                     str(r['end']), r['ref'], r['alt'], \
                                     r['most_severe_impact'], r['aa_change'], \
                                     r['aa_length'], idx_to_sample[idx], \
                                     gts[idx], gene, trans])


def lof_sieve(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

