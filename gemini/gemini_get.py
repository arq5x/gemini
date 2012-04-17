#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import cPickle
import numpy as np
import zlib
#import sql
import sql_extended as sql
from pyparsing import ParseResults
from collections import defaultdict

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
        if not col.startswith("GT_"):
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
    
    # create and execute the main query
    query = "select * from variants"
    if len(main_where) > 0:
        query += " where " + main_where
    # TO DO : handle GROOUP BY, HAVING, ORDER BY, LIMIT


    c.execute(query)
    # (select *) a list of all of the non-gt* columns in the table
    all_cols = [str(tuple[0]) for tuple in c.description if not tuple[0].startswith("gt")]

    # loop through the results of the main query and remove
    # rows that don't meet the conditions placed on the gt_* columns
    gts_requested = False
    if (gts_where != "" or (any("gt_" in s for s in select_cols))):
        gts_requested = True

    for row in c:
        gts       = np.array(cPickle.loads(zlib.decompress(row['gts'])))
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        gt_phases = np.array(cPickle.loads(zlib.decompress(row['gt_phases'])))
        if '*' in select_cols:
            print "\t".join(str(row[col]) for col in all_cols),
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
    def map_samples_to_indicies(c):
        sample_to_idx = {}
        c.execute("select sample_id, name from samples")
        for row in c:
            sample_to_idx[row['name']] = row['sample_id'] - 1
        return sample_to_idx

    sample_to_idx = map_samples_to_indicies(c)
    query_pieces = args.query.split()
    if not any(s.startswith("gt") for s in query_pieces) and \
       not any("gt" in s for s in query_pieces):
       apply_query(c, args.query)
    else:
        (tokens, select_cols, main_where, gts_where) = \
            refine_sql(args.query, sample_to_idx)
        apply_refined_query(c, tokens, select_cols, main_where, gts_where)



def get_shortcut(args, c):
    """
    Router for calling the requested
    shortcut function.
    """
    if args.shortcut == "variants":
        shortcut_variants(args)
    elif args.shortcut == "samples":
        shortcut_samples(args)
    elif args.shortcut == "genotypes":
        shortcut_genotypes(args)
    elif args.shortcut == "tstv":
        shortcut_tstv(args)
    elif args.shortcut == "tstv-coding":
        shortcut_tstv_coding(args)
    elif args.shortcut == "tstv-noncoding":
        shortcut_tstv_noncoding(args)
    elif args.shortcut == "snp-counts":
        shortcut_snpcounts(args)
    elif args.shortcut == "ir-candidates":
        shortcut_ir_candidates(args)
    elif args.shortcut == "sfs":
        shortcut_sfs(args)
    elif args.shortcut == "mds":
        shortcut_mds(c)
    else:
        sys.stderr.write(shortcut + ": unrecognized get shortcut.\n")
        exit()


###############################################################################
# Individual functions for each shortcut.
###############################################################################
def shortcut_genotypes(args):
    query = "SELECT v.chrom, v.start, v.end, \
                     v.ref, v.alt, v.qual, \
                     v.type, v.sub_type, \
                     v.aaf, v.in_dbsnp, v.gene, \
                     v.exonic, v.exon, v.impact, v.is_lof, \
                     s.name, \
                     g.genotype \
             FROM    variants v, \
                     genotypes g, \
                     samples s \
             WHERE   v.variant_id = g.variant_id \
               AND   g.sample_id = s.sample_id"
    sqlite_cli_call(args, query)


def shortcut_ir_candidates(args):
    query = "SELECT v.*, s.name, g.type \
               FROM  variants v, \
                     genotypes g, \
                     samples s \
               WHERE v.variant_id       = g.variant_id \
                AND  g.sample_id       = s.sample_id \
                AND  v.num_hom_alt = 1 \
                AND  v.num_het = 0 \
                AND  g.type   > 0 \
                AND  v.exonic = 1 \
                AND  v.depth >= 200 \
                GROUP BY v.chrom, v.start, v.end, v.impact, v.gene, s.name, g.type" # 1 = het, 2 = hom_alt
    sqlite_cli_call(args, query)


def shortcut_samples(args):
    query = "SELECT * FROM samples"
    sqlite_cli_call(args, query)


def shortcut_snpcounts(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    cmd = "SELECT ref, alt, count(1) \
           FROM   variants \
           WHERE  type = \'snp\' \
           GROUP BY ref, alt"
    # get the ref and alt alleles for all snps.
    c.execute(cmd)

    if args.use_header: args.separator.join(c.description)
    for row in c:
        ref   = str(row[0])
        alt   = str(row[1])
        count = row[2]
        print ref + "->" + alt + args.separator + str(count)


def shortcut_tstv(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    ts_cmd = "SELECT count(1) \
           FROM  variants \
           WHERE type = \'snp\' \
           AND   sub_type = \'ts\'"
    tv_cmd = "SELECT count(1) \
          FROM  variants v \
          WHERE type = \'snp\' \
          AND   sub_type = \'tv\'"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)
    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "transitions" + args.separator + "transversions" + args.separator + "ts/tv"
    print str(ts) + "\t" + str(tv) + "\t" + str(float(ts)/float(tv))


def shortcut_tstv_coding(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    ts_cmd = "SELECT count(1) \
           FROM variants v \
           WHERE v.type = \'snp\' \
           AND v.sub_type = \'ts\' \
           AND v.exonic = 1"
    tv_cmd = "SELECT count(1) \
          FROM variants v \
          WHERE v.type = \'snp\' \
          AND v.sub_type = \'tv\' \
          AND v.exonic = 1"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)
    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "transitions" + args.separator + "transversions" + args.separator + "ts/tv"
    print str(ts) + "\t" + str(tv) + "\t" + str(float(ts)/float(tv))


def shortcut_tstv_noncoding(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    ts_cmd = "SELECT count(1) \
           FROM variants v \
           WHERE v.type = \'snp\' \
           AND v.sub_type = \'ts\' \
           AND v.exonic = 0"
    tv_cmd = "SELECT count(1) \
          FROM variants v \
          WHERE v.type = \'snp\' \
          AND v.sub_type = \'tv\' \
          AND v.exonic = 0"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)

    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "transitions" + args.separator + "transversions" + args.separator + "ts/tv"
    print str(ts) + args.separator + str(tv) + args.separator + str(float(ts)/float(tv))


def shortcut_variants(args):
    query = "SELECT * FROM variants"
    sqlite_cli_call(args, query)


def shortcut_sfs(args):
    query = "SELECT round(aaf," + str(args.precision) + "), count(1) \
             FROM (select aaf from variants group by variant_id) \
             GROUP BY round(aaf," + str(args.precision) +")"
    sqlite_cli_call(args, query)


def shortcut_mds(c):
    """
    Compute the pairwise genetic distance between each sample. 
    """
    idx_to_sample = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        idx_to_sample[int(row['sample_id']) - 1] = row['name']
    
    mds_cmd = "SELECT DISTINCT v.variant_id, v.gt_types\
               FROM variants v\
               WHERE v.type = 'snp'"
    c.execute(mds_cmd)

    # keep a list of numeric genotype values
    # for each sample
    genotypes = defaultdict(list)
    for row in c:
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        for idx, type in enumerate(gt_types):
            genotypes[idx_to_sample[idx]].append(type)
    
    mds = defaultdict(float)
    deno = defaultdict(float)
    # convert the genotype list for each sample
    # to a numpy array for performance.
    # masks stores an array of T/F indicating which genotypes are
    # known (True, [0,1,2]) and unknown (False [-1]). 
    masks = {}
    for sample in genotypes:
        x = np.array(genotypes[sample])
        genotypes[sample] = x
        masks[sample] = \
            np.ma.masked_where(genotypes[sample]>=0, genotypes[sample]).mask
    # compute the euclidean distance for each s1/s2 combination
    # using numpy's vectorized sum() and square() operations.
    # we use the mask arrays to identify the indices of known genotypes
    # for each sample.  by doing a bitwise AND of the mask arrays for the
    # two samples, we have a mask array of variants where __both__ samples
    # were called.
    for s1 in genotypes:
        for s2 in genotypes:
            pair = (s1,s2)
            # which variants have known genotypes for both samples?
            both_mask = masks[s1] & masks[s2]
            gt1 = genotypes[s1]
            gt2 = genotypes[s2]
            eucl_dist = float(np.sum(np.square((gt1-gt2)[both_mask]))) \
                        / \
                        float(np.sum(both_mask))
            mds[pair] = eucl_dist
            deno[pair] = np.sum(both_mask)

    for pair in mds:
        print "\t".join([str(pair), str(mds[pair]/deno[pair])])


def get(parser, args):

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
        elif args.shortcut is not None:
            get_shortcut(args, c)
    else:
        pass

if __name__ == "__main__":
    main()
