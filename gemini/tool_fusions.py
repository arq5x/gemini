#!/usr/bin/env python
import os

import GeminiQuery
from gemini_constants import *
import gemini_subjects as subjects

def report_fusion(event, subjects_dict, args):
    """
    Report the fusion event.
    """
    # filter single line events
    if len(event) == 1:
        sv = event.pop()
        gene1 = sv['gene']
        gene1_strand = sv['strand']
        gene1_start = sv['transcript_min_start']
        gene1_end = sv['transcript_max_end']

        # query the table to test whether the END breakpoint lies in a gene
        gq = GeminiQuery.GeminiQuery(args.db)

        query = """SELECT gene,
                          strand,
                          in_cosmic_census
                   FROM   gene_summary
                   WHERE  gene_summary.chrom = '%s'
                   AND    (gene_summary.transcript_min_start > %s
                          OR gene_summary.transcript_max_end < %s)
                   AND    gene_summary.transcript_min_start < %s
                   AND    gene_summary.transcript_max_end > %s
                   AND    gene_summary.gene != 'None'
                   LIMIT  1
                """ % (sv['chrom'],
                       sv['transcript_max_end'],
                       sv['transcript_min_start'],
                       sv['end'],
                       sv['end'])

        gq.run(query)
        gene2, gene2_strand, gene2_cosmic = (None, None, None)
        for row in gq:
            gene2 = row['gene']
            gene2_strand = row['strand']
            gene2_cosmic = row['in_cosmic_census']
            break # just get the first gene interrupted by the breakend

        # Break if breakpoint2 is intergenic
        if gene2 == None:
            return

        # if SV is a deletion or duplication, genes must be same strand for fusion
        if sv['sub_type'] == 'DEL' or sv['sub_type'] == 'DUP':
            if gene1_strand != gene2_strand:
                return
        # if SV is an inversion, genes must be opposite strands for fusion
        if sv['sub_type'] == 'INV':
            if gene1_strand == gene2_strand:
                return
        # check COSMIC status, if required
        if args.in_cosmic_census and not (sv['in_cosmic_census'] or gene2_cosmic):
            return

        # pass the variables for compatibility with multi-line variants
        end1 = sv
        end2_chrom = end1['chrom']
        end2_start = sv['sv_cipos_start_right']
        end2_end = sv['sv_cipos_end_right']


    # filter multi-line events
    elif len(event) == 2:
        end1 = event.pop()
        end2 = event.pop()
        gene1_strand, gene2_strand = end1['strand'], end2['strand'] # this is gene_summary.strand
        # require that the genes are non-overlapping
        if (end1['chrom'] == end2['chrom'] \
                and end1['transcript_max_end'] >= end2['transcript_min_start'] \
                and end1['transcript_min_start'] <= end2['transcript_max_end']):
            return
        # if breakpoint joins same strand,
        # then genes must be same strand for fusion
        if (end1['sv_strand'][0] == end1['sv_strand'][1] \
                and gene1_strand != gene2_strand):
            return
        # if breakpoint joins opposite strands,
        # then genes must also be opposite strands for fusion
        if (end1['sv_strand'][0] != end1['sv_strand'][1] \
                  and gene1_strand == gene2_strand):
            return
        # check COSMIC status, if required
        if args.in_cosmic_census and not (end1['in_cosmic_census'] or end2['in_cosmic_census']):
            return

        # store the second end for compatibility with single-line SVs
        gene2 = end2['gene']
        end2_chrom = end2['chrom']
        end2_start = end2['sv_cipos_start_right']
        end2_end = end2['sv_cipos_end_right']

    # fusion passes all filters, print
    print '\t'.join(map(str,
                         [end1['chrom'],
                          end1['sv_cipos_start_left'] - 1,
                          end1['sv_cipos_end_left'],
                          end2_chrom,
                          end2_start - 1,
                          end2_end,
                          end1['sv_event_id'],
                          end1['qual'],
                          end1['sv_strand'][0],
                          end1['sv_strand'][1],
                          end1['sub_type'],
                          end1['gene'],
                          gene2,
                          end1['sv_tool'],
                          end1['sv_evidence_type'],
                          end1['sv_is_precise'],
                          ','.join(end1['variant_samples'])
                          ])
                     )
    return

def get_fusions(args):
    """
    Identify candidate rearrangments resulting in fusion genes.
    """
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    idx_to_sample = gq.idx_to_sample
    subjects_dict = subjects.get_subjects(args)

    # create strings for gemini query of command line args
    qual_string, ev_type_string, cosmic_string = ("", "", "")
    if args.min_qual:
        qual_string = " AND qual >= %s" % args.min_qual
    if args.evidence_type:
        ev_type_string = " AND sv_evidence_type = '%s'" % args.evidence_type

    query = """SELECT variants.chrom, start, end,
                      ref, alt,
                      qual,
                      is_somatic, somatic_score,
                      type, sub_type, variants.gene, 
                      sv_strand, sv_length,
                      sv_cipos_start_left,
                      sv_cipos_start_right,
                      sv_cipos_end_left,
                      sv_cipos_end_right,
                      sv_event_id, sv_mate_id,
                      sv_tool, sv_evidence_type,
                      sv_is_precise,
                      gene_summary.strand,
                      gene_summary.transcript_min_start,
                      gene_summary.transcript_max_end,
                      gene_summary.in_cosmic_census
               FROM variants, gene_summary
               WHERE is_somatic = 1   
               AND   type = 'sv'
               AND   variants.gene is not NULL
               AND   variants.chrom = gene_summary.chrom
               AND   variants.gene = gene_summary.gene
               %s
               %s
               ORDER BY sv_event_id
            """ % (qual_string, ev_type_string)

    curr = None
    prev = None
    gq.run(query)
    for row in gq:
        # single-line variants (DEL, DUP, INV)
        if row['sub_type'] != 'complex':
            report_fusion([row], subjects_dict, args)

        # multi-line variants (BND)
        elif row['sv_mate_id']:
            curr = row
            # the SV event ids match, and prev is not None
            if (prev and curr['sv_event_id'] == prev['sv_event_id']):
                report_fusion([prev, curr], subjects_dict, args)
            # shift the previous
            prev = curr

def run(parser, args):
    if os.path.exists(args.db):
        get_fusions(args)
