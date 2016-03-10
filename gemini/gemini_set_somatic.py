#!/usr/bin/env python
from gemini_constants import *
import gemini_subjects
import GeminiQuery

def tag_somatic_mutations(args):

    t_n_pairs = gemini_subjects.get_families(args.db)

    gq = GeminiQuery.GeminiQuery(args.db)

    depth_string, qual_string, ssc_string, chrom_string = ("", "", "", "")
    if args.min_depth:
        depth_string = " AND depth >= %s" % args.min_depth
    if args.min_qual:
        qual_string = " AND qual >= %s" % args.min_qual
    if args.min_somatic_score:
        ssc_string = " AND (type='sv' \
                         OR somatic_score >= %s)" % args.min_somatic_score
    if args.chrom:
        chrom_string = " AND chrom = '%s'" % args.chrom

    if args.chrom is None:
        query = "SELECT variant_id, chrom, start, end, \
                        ref, alt, gene, impact, gts, gt_types, \
                        gt_ref_depths, gt_alt_depths \
                 FROM variants \
                 WHERE 1 \
                 %s \
                 %s \
                 %s \
                 %s" % (depth_string, qual_string, ssc_string, chrom_string)

    gq.run(query)
    smp2idx = gq.sample_to_idx

    somatic_counter = 0
    somatic_v_ids = []

    if args.dry_run:
        print'\t'.join(['tum_name', 'tum_gt', 'tum_alt_freq', 'tum_alt_depth', 'tum_depth', \
                        'nrm_name', 'nrm_gt', 'nrm_alt_freq', 'nrm_alt_depth', 'nrm_depth',
                        'chrom', 'start', 'end', 'ref', 'alt', 'gene'])

    for row in gq:
        # we can skip variants where all genotypes are identical
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

               try:
                   tum_alt_freq = float(tum_alt_depth) / \
                                  (float(tum_alt_depth) + float(tum_ref_depth))
               except ZeroDivisionError:
                   tum_alt_freq = 'NA'

               try:
                   nrm_alt_freq = float(nrm_alt_depth) / \
                                  (float(nrm_alt_depth) + float(nrm_ref_depth))
               except ZeroDivisionError:
                   nrm_alt_freq = 'NA'

               # apply evidence thresholds.
               if (args.max_norm_alt_freq and nrm_alt_freq > args.max_norm_alt_freq) \
                  or \
                  (args.max_norm_alt_count and nrm_alt_depth > args.max_norm_alt_count):
                  continue

               somatic_counter += 1
               somatic_v_ids.append((1, row['variant_id']))

               print'\t'.join(str(s) for s in [tumor.name,  tum_gt, tum_alt_freq, tum_alt_depth, tum_depth, \
                                   normal.name, nrm_gt, nrm_alt_freq, nrm_alt_depth, nrm_depth, \
                                   row['chrom'], row['start'], row['end'], row['ref'], row['alt'], row['gene']])

    if not args.dry_run:
        import database
        conn, metadata = database.get_session_metadata(args.db)

        # now set the identified mutations to True.
        update_qry = "UPDATE variants SET is_somatic = 1 "
        update_qry += " WHERE variant_id IN (%s)"
        update_qry %= ",".join(str(x[1]) for x in somatic_v_ids)
        res = conn.execute(update_qry)
        assert res.rowcount == somatic_counter
        print "Identified and set", somatic_counter, "somatic mutations"
        conn.commit()
    else:
        print "Would have identified and set", somatic_counter, "somatic mutations"

def set_somatic(parser, args):

    tag_somatic_mutations(args)
