#!/usr/bin/env python
import sqlite3
import numpy as np
import zlib
import re
import os
import cPickle
import gemini_utils as util

import pybedtools as pbt

def get_window_data(c, analysis_type):
    """
    Create a temp file of the requested statistic for each variant.
    
    Execute a query against the variants table
    that extracts the requested column for each variant.
    save the results to '.temp', which will be loaded
    into a pybedtools BedTool for use with the bedtools map 
    function.  This will compute the requested statistic
    for each variant in the variants table
    """
    if analysis_type == "hwe": column = 'hwe'
    elif analysis_type == "nucl_div": column = 'pi'
    
    t = open('.temp', 'w')
    query = "SELECT chrom,start,end," + \
             column + \
             " FROM variants ORDER BY chrom,start"
    c.execute(query)
    for row in c:
        t.write('%s\t%d\t%d\t%f\n' % (str(row['chrom']),
                                   int(row['start']),
                                   int(row['end']),
                                   float(row[column])))
    t.close()
    # Tell bedtools map that the statistuc is in the fourth column.
    # Parameterized for future mods,
    return 4


def make_windows(c, args):
    """
    Compute the requested statistic for the user-defined windows.
    """    
    # create our windows with pybedtools
    window = pbt.BedTool()
    windows = window.window_maker(genome='hg19', 
                                  w=args.window_size, 
                                  s=args.step_size)
    
    # create a temp file ('.temp') storing the requested stat 
    # for each variant. Load this into a pybedtools BedTool
    op_col = get_window_data(c, args.analysis_type)
    window_data = pbt.BedTool('.temp')

    # Use bedtools map to summarize and report
    # the requested statistic for each window
    windowed_analysis = windows.map(window_data, o=args.op_type, c=op_col)
    for window in windowed_analysis:
        print window,
    # cleanup
    os.remove('.temp')
       
def windower(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        # on y va
        make_windows(c, args)
