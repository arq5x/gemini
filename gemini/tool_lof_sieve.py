#!/usr/bin/env python
from __future__ import absolute_import, print_function

import re
import sqlalchemy as sql
from . import gemini_utils as util
from .gemini_constants import *
from . import compression as Z


def get_ind_lof(conn, metadata, args):

    idx_to_sample = util.map_indices_to_samples(metadata)

    query = "SELECT v.chrom, v.start, v.end, v.ref, v.alt, \
                             v.impact, v.aa_change, v.aa_length, \
                             v.gt_types, v.gts, i.gene, \
                             i.transcript,  i.biotype\
             FROM variants v, variant_impacts i \
             WHERE v.variant_id = i.variant_id \
             AND i.is_lof='1' \
             AND v.type = 'snp'"

    res = conn.execute(sql.text(query))

    # header
    print('\t'.join(['chrom', 'start', 'end', 'ref', 'alt',
                     'highest_impact', 'aa_change', 'var_trans_pos',
                     'trans_aa_length', 'var_trans_pct',
                     'sample', 'genotype', 'gene', 'transcript',
                     'trans_type']))

    unpack = Z.unpack_genotype_blob
    for r in res:
        try:
            gt_types = unpack(r['gt_types'])
            gts = unpack(r['gts'])
        except:
            unpack = Z.snappy_unpack_blob
            gt_types = unpack(r['gt_types'])
            gts = unpack(r['gts'])

        gene = str(r['gene'])
        trans = str(r['transcript'])

        aa_change = str(r['aa_change'])
        aa_length = str(r['aa_length'])
        transcript_pos = None
        transcript_pct = None
        if aa_change != 'None':
            try:
                #transcript_pos for snpEff annotated VCF
                transcript_pos = re.findall('\S(\d+)\S', aa_change)[0]
            except IndexError:
                #transcript_pos for VEP annotated VCF
                if aa_length != 'None' and \
                        aa_length.split("/")[0] != "-":
                    transcript_pos = aa_length.split("/")[0]
        #handle non exonic variants
        if transcript_pos is None:
            transcript_pct = '/'
        #transcript_pct for snpEff annotated VCF
        elif aa_length != 'None' and "/" not in aa_length:
            transcript_pct = float(transcript_pos) / float(aa_length)
        #transcript_pct for VEP annotated VCF
        elif aa_length != 'None' and "/" in aa_length:
            transcript_pct = float(transcript_pos) / float(aa_length.split("/")[1])

        for idx, gt_type in enumerate(gt_types):
            if gt_type == HET or gt_type == HOM_ALT:
                print("\t".join([r['chrom'], str(r['start']),
                                 str(r['end']), r['ref'], r['alt'],
                                 r['impact'],
                                 r['aa_change'] or 'None',
                                 transcript_pos or 'None',
                                 r['aa_length'] or 'None',
                                 str(transcript_pct) or 'None',
                                 idx_to_sample[idx],
                                 gts[idx], gene, trans, r['biotype'] or
                                 'None']))


def lof_sieve(parser, args):

    from . import database
    conn, metadata = database.get_session_metadata(args.db)
    get_ind_lof(conn, metadata, args)
