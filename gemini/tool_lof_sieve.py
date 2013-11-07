#!/usr/bin/env python

import os
import re
import sqlite3
import gemini_utils as util
from gemini_constants import *


def get_ind_lof(c, args):

    idx_to_sample = util.map_indices_to_samples(c)

    query = "SELECT v.chrom, v.start, v.end, v.ref, v.alt, \
                             v.impact, v.aa_change, v.aa_length, \
                             v.gt_types, v.gts, i.gene, \
                             i.transcript,  i.biotype\
             FROM variants v, variant_impacts i \
             WHERE v.variant_id = i.variant_id \
             AND i.is_lof='1' \
             AND v.type = 'snp'"

    c.execute(query)

    # header
    print '\t'.join(['chrom', 'start', 'end', 'ref', 'alt',
                     'highest_impact', 'aa_change', 'var_trans_pos',
                     'trans_aa_length', 'var_trans_pct',
                     'sample', 'genotype', 'gene', 'transcript', 'trans_type'])

    for r in c:
        gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
        gts = np.array(cPickle.loads(zlib.decompress(r['gts'])))
        gene = str(r['gene'])
        trans = str(r['transcript'])

        aa_change = str(r['aa_change'])
        aa_length = str(r['aa_length'])
        transcript_pos = None
        transcript_pct = None
        if aa_change != 'None':
            transcript_pos = re.findall('\S(\d+)\S', aa_change)[0]
            if aa_length != 'None':
                transcript_pct = float(transcript_pos) / float(aa_length)

        for idx, gt_type in enumerate(gt_types):
            if gt_type == HET or gt_type == HOM_ALT:
                print "\t".join([r['chrom'], str(r['start']),
                                 str(r['end']), r['ref'], r['alt'],
                                 r['impact'],
                                 r['aa_change'] or 'None',
                                 transcript_pos or 'None',
                                 r['aa_length'] or 'None',
                                 str(transcript_pct) or 'None',
                                 idx_to_sample[idx],
                                 gts[idx], gene, trans, r['biotype']])


def lof_sieve(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        get_ind_lof(c, args)
