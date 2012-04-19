#!/usr/bin/env python
import sqlite3
import os
import sys
import numpy as np
import cPickle
import zlib
import collections

import gemini_utils as util
from nesteddict import NestedDict

def get_compound_hets(c, args):
    """
    Report the count of each genotype class
    observed for each sample.
    """
    idx_to_sample = util.map_indicies_to_samples(c)

    comp_hets = NestedDict()

    query = "SELECT * FROM variants where is_exonic = 1"
    c.execute(query)
    
    # count the number of each genotype type obs. for each sample.
    for row in c:
        gene      = row['gene']
        chrom     = row['chrom']
        start     = row['start']
        end       = row['end']
        ref       = row['ref']
        alt       = row['alt']
        impact    = row['impact']
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        for idx, gt_type in enumerate(gt_types):
            if gt_type == 1: 
                comp_hets[sample][gene].append(start)
    # report
    for sample in comp_hets:
        for gene in comp_hets[sample]:
            print sample, gene, comp_hets[sample][gene]


def run(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        
        get_compound_hets(c, args)

