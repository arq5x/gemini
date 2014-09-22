#!/usr/bin/env python
import sqlite3
import os
import sys
from collections import defaultdict

import GeminiQuery
import gemini_utils as util
from gemini_constants import *
import gemini_subjects as subjects

def report_fusion(event):
    """
    Report the fusion event.
    """
    # Colby, this is where you could apply any further exclusion criteria, etc.
    end1 = event.pop()
    end2 = event.pop()
    
    # if exclusion:
    #   continue
    # if exclusion:
    #   continue
    print "--->fusion gene!"
    print end1, "\n", end2

def get_fusions(args):
    """
    Identify candidate rearrangments resulting in fusion genes.
    """
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    idx_to_sample = gq.idx_to_sample
    subjects_dict = subjects.get_subjects(args)

    ##### COLBY, change "WHERE is_somatic is NULL"  to "WHERE is_somatic = 1" 
    query = """SELECT chrom, start, end, 
                      is_somatic, somatic_score,
                      type, sub_type, gene, 
                      sv_strand, sv_length,
                      sv_event_id, sv_mate_id,
                      sv_tool, sv_evidence_type
               FROM variants
               WHERE is_somatic is NULL    
               AND   sv_mate_id is not NULL
               AND   type = 'sv'
               AND   sub_type = 'complex'
               AND   gene is not NULL
               ORDER BY sv_event_id
            """

    curr = None
    prev = None
    events = []
    gq.run(query)
    for row in gq:
        curr = row['sv_event_id']
        # the SV event id changed.
        if curr != prev and prev is not None:
            # did both ends of the sv meet all the query criteria
            # and are both ends on the same strand?
            if len(events) == 2 and \
                (events[0]['sv_strand'] == events[1]['sv_strand']):
                report_fusion(events)
            # we are done with this candidate
            events = []   
        else:
            events.append(row)
        prev = curr

def run(parser, args):
    if os.path.exists(args.db):
        get_fusions(args)
        
