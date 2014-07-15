#!/usr/bin/env python
import sqlite3
from collections import Counter
from gemini_constants import *
import gemini_subjects
import GeminiQuery

def tag_somatic_mutations(args):

    t_n_pairs = gemini_subjects.get_families(args.db)

    gq = GeminiQuery.GeminiQuery(args.db)

    if args.chrom is None:
        query = "SELECT variant_id, chrom, start, end, \
                        ref, alt, gene, impact, gts, gt_types, \
                        gt_ref_depths, gt_alt_depths \
                 FROM variants \
                 WHERE depth >= " + str(args.min_depth) + \
                 " AND   qual >= " + str(args.min_qual) 
    else:
        query = "SELECT variant_id, chrom, start, end, \
                ref, alt, gene, impact, gts, gt_types, \
                gt_ref_depths, gt_alt_depths \
         FROM variants \
         WHERE depth >= " + str(args.min_depth) + \
         " AND qual >= " + str(args.min_qual) + \
         " AND chrom = \'" + args.chrom + "\'"

    gq.run(query)
    smp2idx = gq.sample_to_idx

    somatic_counter = 0
    somatic_v_ids = []

    if args.dry_run:
        print'\t'.join(['tum_name', 'tum_gt', 'tum_alt_freq', 'tum_alt_depth', 'tum_depth', \
                        'nrm_name', 'nrm_gt', 'nrm_alt_freq', 'nrm_alt_depth', 'nrm_depth',
                        'chrom', 'start', 'end', 'ref', 'alt', 'gene'])

    for row in gq:

        # we can skip varinats where all genotypes are identical
        if len(set(row['gt_types'])) == 1:
            continue

        for pair in t_n_pairs:

            samples = pair.subjects
            if len(samples) != 2:
                continue

            tumor = pair.subjects[0]
            normal = pair.subjects[1]
            # swap if we guessed the tumor incorrectly
            if tumor.affected is False:
                tumor, normal = normal, tumor

            tum_idx = smp2idx[tumor.name]
            nrm_idx = smp2idx[normal.name]

            tum_gt = row['gts'][tum_idx]
            nrm_gt = row['gts'][nrm_idx]

            tum_gt_type = row['gt_types'][tum_idx]
            nrm_gt_type = row['gt_types'][nrm_idx]

            if nrm_gt_type == tum_gt_type:
                continue

            if nrm_gt_type == UNKNOWN or tum_gt_type == UNKNOWN:
                continue

            # the genotypes pass the smell test for somatic
            # mutations if in this block.
            if (nrm_gt_type == HOM_REF and tum_gt_type != HOM_REF):

               
               tum_ref_depth = row['gt_ref_depths'][tum_idx]
               nrm_ref_depth = row['gt_ref_depths'][nrm_idx]

               tum_alt_depth = row['gt_alt_depths'][tum_idx]
               nrm_alt_depth = row['gt_alt_depths'][nrm_idx]

               # total observed depth
               nrm_depth = nrm_alt_depth + nrm_ref_depth
               tum_depth = tum_alt_depth + tum_ref_depth

               if (nrm_depth < args.min_norm_depth \
                  or \
                  tum_depth < args.min_tumor_depth):
                  continue

               tum_alt_freq = float(tum_alt_depth) / \
                              (float(tum_alt_depth) + float(tum_ref_depth))

               nrm_alt_freq = float(nrm_alt_depth) / \
                              (float(nrm_alt_depth) + float(nrm_ref_depth))

               # apply evidence thresholds.
               if nrm_alt_freq > args.max_norm_alt_freq \
                  or \
                  nrm_alt_depth > args.max_norm_alt_count:
                  continue

               somatic_counter += 1
               somatic_v_ids.append((1, row['variant_id']))
               
               print'\t'.join(str(s) for s in [tumor.name,  tum_gt, tum_alt_freq, tum_alt_depth, tum_depth, \
                                   normal.name, nrm_gt, nrm_alt_freq, nrm_alt_depth, nrm_depth, \
                                   row['chrom'], row['start'], row['end'], row['ref'], row['alt'], row['gene']])

    if not args.dry_run:
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        c = conn.cursor()

        # now set the identified mutations to True.
        update_qry = "UPDATE variants SET is_somatic = ? "
        update_qry += " WHERE variant_id = ?"
        c.executemany(update_qry, somatic_v_ids)
        print "Identified and set", somatic_counter, "somatic mutations"
    else:
        print "Would have identified and set", somatic_counter, "somatic mutations"

def set_somatic(parser, args):
    
    tag_somatic_mutations(args)

