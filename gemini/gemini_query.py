#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import cPickle
import numpy as np
import zlib
import string
import itertools
from pyparsing import ParseResults
from collections import defaultdict
from collections import OrderedDict

# gemini imports
import gemini_utils as util
from gemini_constants import *
from gemini_utils import OrderedSet as oset
import sql_extended as sql
import compression

# REGEX to trap the special gt_* columns
gt_sel_re = re.compile("gt(\S*)\.(\S+)")

def _flatten(l):
    return [item for sublist in l for item in sublist]

def _correct_genotype_col(raw_col, sample_to_idx):
    """
    Convert a _named_ genotype index to a _numerical_
    genotype index so that the appropriate value can be
    extracted for the sample from the genotype numpy arrays.
    
    These lookups will be eval()'ed on the resuting rows to
    extract the appropriate information.
    
    For example, convert gt_types.1478PC0011 to gt_types[11]
    """
    if raw_col == "*":
       return raw_col.lower()
    (column, sample) = raw_col.split('.')
    corrected = column.lower() + "[" + str(sample_to_idx[sample]).lower() + "]"
    return corrected


def _split_select(query, sample_to_idx):
    """
    Build a list of _all_ columns in the SELECT statement
    and segregated the non-genotype specific SELECT columns.
    
    This is used to control how to report the results, as the
    genotype-specific columns need to be eval()'ed whereas others
    do not.
    
    For example: "SELECT chrom, start, end, gt_types.1478PC0011"
    will populate the lists as follows:
    
    select_columns = ['chrom', 'start', 'end']
    all_columns = ['chrom', 'start', 'end', 'gt_types[11]']
    """
    select_columns = []
    all_columns_new = []
    all_columns_orig = []
    gt_col_map = {}

    # iterate through all of the select columns andclear
    # distinguish the genotype-specific columns from the base columns
    from_loc = query.find("from")
    if from_loc < 1:
         sys.exit("Malformed query.")

    raw_select_clause = query[0:from_loc].rstrip()
    rest_of_query = query[from_loc:len(query)]

    for token in raw_select_clause.replace(',','').split():
        if token == "SELECT" or token == "select":
            continue
        if not token.startswith("GT") and not token.startswith("gt"):
            select_columns.append(token)
            all_columns_new.append(token)
            all_columns_orig.append(token)
        else:
            new_col = _correct_genotype_col(token, sample_to_idx)
            all_columns_new.append(new_col)
            all_columns_orig.append(token)
            gt_col_map[token] = new_col

    return select_columns, all_columns_new, all_columns_orig, gt_col_map


def add_gt_cols_to_query(query):
    """
    We have to modify the raw query to select the genotype
    columns in order to support the genotype filters.  That is,
    if the user wants to limit the rows returned based upon, for example,
    "gts.joe == 1", then we need to select the full gts BLOB column in
    order to enforce that limit.  The user wouldn't have selected gts as a
    columns, so therefore, we have to modify the select statement to add
    it.
    
    In essence, when a gneotype filter has been requested, we always add
    the gts, gt_types and gt_phases columns.
    """
    from_loc = query.find("from")
    if from_loc > 1:
        raw_select_clause = query[0:from_loc].rstrip()
        rest_of_query = query[from_loc:len(query)]
        
        # remove any GT columns
        select_clause_list = []
        for token in raw_select_clause.split():
            if not token.startswith("gt") and not token.startswith("GT"):
                select_clause_list.append(token)

        # add the genotype columns to the query
        if select_clause_list[len(select_clause_list)-1].endswith(",") or \
               (len(select_clause_list) == 1 and \
               select_clause_list[0].strip().lower() == "select"):
            select_clause = " ".join(select_clause_list) + \
                                 " gts, gt_types, gt_phases, gt_depths "
        else:
            select_clause = " ".join(select_clause_list) + \
                                 ", gts, gt_types, gt_phases, gt_depths "
        query =select_clause + rest_of_query
        # extract the original select columns
        return query
        
    else:
       sys.exit("Malformed query.")


def apply_basic_query(c, query, use_header):
    """
    Execute a vanilla query. That is, not gt* columns
    are in either the select or where clauses.
    """
    c.execute(query)
    all_cols = [str(tuple[0]) for tuple in c.description \
                                    if not tuple[0].startswith("gt")]

    if use_header:
        h = [col for col in all_cols]
        yield OrderedDict(itertools.izip(h,h))

    for row in c:
        yield OrderedDict(row)


def apply_query_w_genotype_select(c, query, use_header):
    """
    Execute a query that contains gt* columns in only in the SELECT.
    """
    # construct a mapping of sample names to list indices
    sample_to_idx = util.map_samples_to_indicies(c)
    
    (select_cols, all_cols_new, all_cols_orig, gt_col_map) = \
                                    _split_select(query, sample_to_idx)
    
    query = add_gt_cols_to_query(query.lower())
    c.execute(query)
    
    # what are the columns that were actually selected by the user.
    all_query_cols = [str(tuple[0]) for tuple in c.description \
                              if not tuple[0].startswith("gt")]

    if "*" in select_cols:
        select_cols.remove("*")
        #all_cols_orig.remove("*")
        all_cols_new.remove("*")
        select_cols += all_query_cols
        
    if use_header:
        h = [col for col in all_query_cols] + \
            [col for col in oset(all_cols_orig) - oset(select_cols)]
        yield OrderedDict(itertools.izip(h,h))

    report_cols = all_query_cols + list(oset(all_cols_new) - oset(select_cols))
    for row in c:
        gts       = compression.unpack_genotype_blob(row['gts'])
        gt_types  = compression.unpack_genotype_blob(row['gt_types'])
        gt_phases = compression.unpack_genotype_blob(row['gt_phases'])
        gt_depths = compression.unpack_genotype_blob(row['gt_depths'])

        fields = OrderedDict()
        for idx, col in enumerate(report_cols):
            if col == "*": continue
            if not col.startswith("gt") and not col.startswith("GT"):
                fields[col] = row[col]
            else:
                fields[col] = eval(col.strip())
        yield fields


def get_query_file(args):
    """
    Execute a batch of user-defined queries
    passed in as a file from the command line.
    """
    query = ""
    for line in open(query_file, 'r'):
        query += line.strip() + " "
        sqlite_cli_call(args, query)


def get_query(args, c):
    """
    Execute a user-defined query passed in via
    the command line.
    """
    sample_to_idx = util.map_samples_to_indicies(c)
    
    query_pieces = args.query.split()
    if not any(s.startswith("gt") for s in query_pieces) and \
       not any("gt" in s for s in query_pieces):
        apply_basic_query(c, args)
    else:
        apply_query_w_genotype_select(c, args.query, args.use_header)


def filter_query(c, query, gt_filter, use_header):
    """
    Execute a base SQL query while applying filters on the returned 
    rows based on filters applied to the genotype-specific columns.
    
    For example:
    --gt_filter "(gt_types.1478PC0011 == 1 or gt_types.1478PC0012 == 1)
    """

    def correct_genotype_filter(gt_filter, sample_to_idx):
        """
        This converts a "raw" genotype filter supplied by the user
        to a filter than can be eval()'ed.  Specifically, we must
        convery a _named_ genotype index to a _numerical_
        genotype index so that the appropriate value can be
        extracted for the sample from the genotype numpy arrays.
        
        For example, converts:
        --gt-filter "(gt_types.1478PC0011 == 1)"
        to
        (gt_types[11] == 1)
        """
        corrected_gt_filter = []
        tokens = re.split(r'[\s+]+', gt_filter)
        for token in tokens:
            if token.find("gt") >= 0 or token.find("GT") >= 0:
                corrected = _correct_genotype_col(token, sample_to_idx)
                corrected_gt_filter.append(corrected)
            else:
                corrected_gt_filter.append(token)
        return " ".join(corrected_gt_filter)


    # construct a mapping of sample names to list indices
    sample_to_idx = util.map_samples_to_indicies(c)
    
    gt_filter = correct_genotype_filter(gt_filter, sample_to_idx)
    (select_cols, all_cols_new, all_cols_orig, gt_col_map) = \
                                    _split_select(query, sample_to_idx)

    query = add_gt_cols_to_query(query.lower())
    
    c.execute(query)
    
    # what are the columns that were actually selected by the user.
    all_query_cols = [str(tuple[0]) for tuple in c.description \
                              if not tuple[0].startswith("gt")]

    if "*" in select_cols:
        select_cols.remove("*")
        all_cols_orig.remove("*")
        all_cols_new.remove("*")
        select_cols += all_query_cols
    
    if use_header:
        h = [col for col in all_query_cols] + \
            [col for col in oset(all_cols_orig) - oset(select_cols)]
        yield OrderedDict(itertools.izip(h,h))


    report_cols = all_query_cols + list(oset(all_cols_new) - oset(select_cols))
    for row in c:
        gts       = compression.unpack_genotype_blob(row['gts'])
        gt_types  = compression.unpack_genotype_blob(row['gt_types'])
        gt_phases = compression.unpack_genotype_blob(row['gt_phases'])
        gt_depths = compression.unpack_genotype_blob(row['gt_depths'])

        if not eval(gt_filter):
            continue
        
        fields = OrderedDict()
        for idx, col in enumerate(report_cols):
            if col == "*": continue
            if not col.startswith("gt") and not col.startswith("GT"):
                fields[col] = row[col]
            else:
                fields[col] = eval(col.strip())
        yield fields


def query(parser, args):

    if (args.db is None):
        parser.print_help()
    # open up a new database
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row # allow us to refer to columns by name
        c = conn.cursor()
        
        row_iter = None
        if args.query is not None:
            if args.gt_filter is None:
                
                query_pieces = args.query.split()
                if not any(s.startswith("gt") for s in query_pieces) and \
                   not any("gt" in s for s in query_pieces):
                    row_iter = apply_basic_query(c, args.query, args.use_header)
                else:
                    row_iter = apply_query_w_genotype_select(c, \
                                                  args.query, args.use_header)
            else:
                row_iter = filter_query(c, args.query, \
                                        args.gt_filter, args.use_header)
                                    
            for row in row_iter:
                print args.separator.join([str(row[c]) for c in row])
                
        elif args.queryfile is not None:
            get_query_file(args)

if __name__ == "__main__":
    main()
