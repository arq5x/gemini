#!/usr/bin/env python

import os
import sys
import sqlite3
import subprocess
import collections
import pybedtools as pbt
import pysam
import time

import gemini_constants
from gemini.annotations import annotations_in_region, guess_contig_naming

def add_requested_column(col_name, update_cursor):
    """
    Attempt to add a new, user-defined column to the
    variants table.  Exit if the column already exists.
    """
    try:
        alter_qry = "ALTER TABLE variants ADD COLUMN " + col_name + " BOOL DEFAULT NULL"
        update_cursor.execute(alter_qry)
    except sqlite3.OperationalError:
        pass
        #sys.exit("ERROR: column \"" + col_name + "\" already exists in variants table")  

def _annotate_variants(args, conn, get_val_fn):
    """Generalized annotation of variants with a new column.
    
    get_val_fn takes a list of annotations in a region and returns
    the value for that region to update the database with.

    Separates selection and identification of values from update, to avoid concurrent
    database access errors from sqlite3, especially on NFS systems. The retained
    to_update list is small, but batching could help if memory issues emerge.
    """
    # For each, use Tabix to detect overlaps with the user-defined
    # annotation file.  Update the variant row with T/F if overlaps found.
    annos = pysam.Tabixfile(args.anno_file)
    naming = guess_contig_naming(annos)
    select_cursor = conn.cursor()
    select_cursor.execute("SELECT chrom, start, end, variant_id FROM variants")
    to_update = []
    for row in select_cursor:
        to_update.append((str(row["variant_id"]),
                          get_val_fn(annotations_in_region(row, annos, "tuple", naming))))
    update_cursor = conn.cursor()
    add_requested_column(args.col_name, update_cursor)
    for variant_id, val in to_update:
        update_qry = "UPDATE variants SET " + args.col_name + " = " + str(val) + \
                     " WHERE variant_id = " + variant_id
        update_cursor.execute(update_qry)

def annotate_variants_bool(args, conn):
    """
    Populate a new, user-defined column in the variants
    table with a BOOLEAN indicating whether or not
    overlaps were detected between the variant and the 
    annotation file.
    """
    def has_anno_hit(hits):
        has_hit = 0
        for hit in hits:
            has_hit = 1
            break
        return has_hit
    return _annotate_variants(args, conn, has_anno_hit)

def annotate_variants_count(args, conn):
    """
    Populate a new, user-defined column in the variants
    table with a INTEGER indicating the count of overlaps
    between the variant and the 
    annotation file.
    """
    def get_hit_count(hits):
        count = 0
        for hit in hits:
            count += 1
        return count
    return _annotate_variants(args, conn, get_hit_count)

def annotate_variants_list(args, conn):
    """
    Populate a new, user-defined column in the variants
    table with a INTEGER indicating the count of overlaps
    between the variant and the 
    annotation file.
    """
    def get_hit_list(hits):
        hit_list = []
        for hit in hits:
            try:
                hit_list.append(hit[int(args.col_extract) - 1])
            except IndexError:
                sys.exit("Column " + args.col_extract + " exceeds \
                          the number of columns in your \
                          annotation file. Exiting.")
        if len(hit_list) > 0:
            val = ",".join(hit_list)
            return "'%s'" % val
        else:
            return "NULL"
    return _annotate_variants(args, conn, get_hit_list)

def annotate(parser, args):

    if (args.db is None):
        parser.print_help()
        exit()

    if not os.path.exists(args.db):
        sys.stderr.write("Error: cannot find database file.")
        exit(1)
    if not os.path.exists(args.anno_file):
        sys.stderr.write("Error: cannot find annotation file.")
        exit(1)

    conn = sqlite3.connect(args.db)
    conn.row_factory = sqlite3.Row # allow us to refer to columns by name
    conn.isolation_level = None

    if args.col_type == "boolean":
        annotate_variants_bool(args, conn)
    elif args.col_type == "count":
        annotate_variants_count(args, conn) 
    elif args.col_type == "list":
        if args.col_extract is None:
            sys.exit("You must specify which column to extract (-e) from the annotation file.")
        else:
            annotate_variants_list(args, conn)
    else:
        sys.exit("Unknown column type requested. Exiting.")


if __name__ == "__main__":
    main()
