#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import cPickle
import numpy as np
import zlib
from pyparsing import ParseResults
from collections import defaultdict

# gemini imports
import gemini_utils as util
import sql_extended as sql

# REGEX to trap the special gt_* columns
gt_sel_re = re.compile("GT(\S*)\.(\S+)")

def flatten(l):
    return [item for sublist in l for item in sublist]

def refine_sql(query, sample_to_idx):
    def correct_genotype_col(raw_col):
        if raw_col == "*":
           return raw_col.lower()
        match = gt_sel_re.findall(raw_col)[0]
        return "gt" + match[0].lower() + "[" + str(sample_to_idx[match[1]]).lower() + "]"
    
    # used to identify strings that are where condition constructs
    where_keywords = ['where', 'and', 'or']
    
    # use pyparsing in sql.py module to parse raw SQL statement
    tokens = sql.parse_sql(query)
    # adjust any of the gt_* columns to use sample offsets, not names
    select_columns = []
    
    # build a list of the select columns while converting 
    # the genotype columns (GT_*) to their appropriate numpy indices
    for col in tokens.select:
        if not col.startswith("GT"):
            select_columns.append(col)
        else:
            select_columns.append(correct_genotype_col(col))

    # Separate the where clause for the main SQL statement
    # from the conditions placed on the gt* columns.  The
    # latter must be applied post hoc to the rows returned 
    # from the former, as SQL where conditions cannot be applied to BLOBs
    main_where = []
    gtypes_where = []
    last_keyword = None
    is_first_where = True
    # tokens.where[0] b/c pyparsing grammar excessively nests results
    for where_piece in tokens.where:
        # and | or | in
        if where_piece in where_keywords:
            last_keyword = where_piece
            continue

        # at this point, we know we are dealing with an actual condition, so
        # we flatten the pieces of the condition in cases of 
        # nested boolean logic
        # add the condition as normal if not a gt* field
        if not any("GT" in s for s in where_piece):
            if not is_first_where:
                main_where.append(last_keyword)
                is_first_where = False
            for piece in where_piece:
                main_where.append(piece)
        # gt* fileds must be converted from referring to sample name
        # to instead referring to sample offsets in the relevant col. BLOB
        else:
            if not is_first_where:
                gtypes_where.append(last_keyword)
                is_first_where = False
            for piece in where_piece:
                # SQL allows "=" for equality, Python eval needs "=="
                if piece == "=": 
                    gtypes_where.append("==")
                elif piece.startswith("GT"):
                    gtypes_where.append(correct_genotype_col(piece))
                else:
                    gtypes_where.append(piece)

    # hand off the select columns and the two different
    # where conditions.
    return(tokens, select_columns,
           " ".join(main_where),
           " ".join(gtypes_where))


def apply_query(c, query):
    """
    Execute a vanilla query. That is, not gt* columns
    are in either the select or where clauses.
    """
    c.execute(query)
    # (select *) a list of all of the non-gt* columns in the table
    all_cols = [str(tuple[0]) for tuple in c.description if not tuple[0].startswith("gt")]
    for row in c:
        print "\t".join(str(row[col]) for col in all_cols)


def apply_refined_query(c, tokens, select_cols, main_where, gts_where):
    """
    Execute a query that contains gt* columns in either the
    select or the where.
    """
    # construct and execute the main query
    query = "select * from variants"
    if len(main_where) > 0:
        query += " where " + main_where
    # TO DO : handle GROOUP BY, HAVING, ORDER BY, LIMIT

    c.execute(query)
    # (select *) a list of all of the non-gt* columns in the table
    all_cols = [str(tuple[0]) for tuple in c.description if not tuple[0].startswith("gt")]

    # track whether gt_* columns are involved in the SELECT, WHERE or both
    gts_where_req = False
    gts_select_req = False
    if (gts_where != ""):
        gts_where_req = True
    if (any("gt" in s for s in select_cols)):
        gts_select_req = True

    for row in c:
        gts       = np.array(cPickle.loads(zlib.decompress(row['gts'])))
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        gt_phases = np.array(cPickle.loads(zlib.decompress(row['gt_phases'])))
        
        # only report this row if =:
        #    a) there was no where clause for gt_* cols
        # or b) there was a gt_* clause and this row passed
        if gts_where == "" or (gts_where != "" and eval(gts_where)):
            # select *
            if '*' in select_cols:
                print "\t".join(str(row[col]) for col in all_cols),
            # select chrom, is_lof
            elif not gts_select_req:
                print "\t".join(str(row[col]) for col in select_cols),
            # select chrom, gt_types.HG00331, gt_types.HG00332
            else:
                for col in select_cols:
                    if col == "*": continue
                    if not col.startswith("gt"):
                        print str(row[col]) + "\t",
                    else:
                        # e.g., eval gt_types[141] and print result (0,1,2,etc.)
                        print str(eval(col.strip())) + "\t",
            print


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
       apply_query(c, args.query)
    else:
        (tokens, select_cols, main_where, gts_where) = \
            refine_sql(args.query, sample_to_idx)
        apply_refined_query(c, tokens, select_cols, main_where, gts_where)


def query(parser, args):

    if (args.db is None):
        parser.print_help()
    # open up a new database
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row # allow us to refer to columns by name
        c = conn.cursor()
        
        if args.query is not None:
            get_query(args, c)
        elif args.queryfile is not None:
            get_query_file(args)

if __name__ == "__main__":
    main()
