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


def annotate_variants_bool(args, select_cursor, update_cursor):
    """
    Populate a new, user-defined column in the variants
    table with a BOOLEAN indicating whether or not
    overlaps were detected between the variant and the 
    annotation file.
    """
    add_requested_column(args.col_name, update_cursor)

    # For each, use Tabix to detect overlaps with the user-defined
    # annotation file.  Update the variant row with T/F if overlaps found.
    annos = pysam.Tabixfile(args.anno_file)
    naming = guess_contig_naming(annos)
    select_cursor.execute("SELECT chrom, start, end, variant_id FROM variants")
    for row in select_cursor:
        has_hit = False
        for hit in annotations_in_region(row, annos, naming=naming):
            has_hit = True
            break
        if has_hit:
            update_qry = "UPDATE variants SET " + args.col_name + " = 1 " + \
                         " WHERE variant_id = " + str(row['variant_id'])
        else:
            update_qry = "UPDATE variants SET " + args.col_name + " = 0 " + \
                         " WHERE variant_id = " + str(row['variant_id'])
        update_cursor.execute(update_qry)


def annotate_variants_count(args, select_cursor, update_cursor):
    """
    Populate a new, user-defined column in the variants
    table with a INTEGER indicating the count of overlaps
    between the variant and the 
    annotation file.
    """
    add_requested_column(args.col_name, update_cursor)

    # For each, use Tabix to count overlaps with the user-defined
    # annotation file.  Update the variant row with the count.
    annos = pysam.Tabixfile(args.anno_file)
    naming = guess_contig_naming(annos)
    select_cursor.execute("SELECT chrom, start, end, variant_id FROM variants")
    for row in select_cursor:
        count = 0
        for hit in annotations_in_region(row, annos, naming=naming):
            count += 1
        update_qry = "UPDATE variants SET " + args.col_name + " = " + str(count) + \
                     " WHERE variant_id = " + str(row['variant_id'])
        update_cursor.execute(update_qry)

def annotate_variants_list(args, select_cursor, update_cursor):
    """
    Populate a new, user-defined column in the variants
    table with a INTEGER indicating the count of overlaps
    between the variant and the 
    annotation file.
    """
    add_requested_column(args.col_name, update_cursor)

    # For each, use Tabix to count overlaps with the user-defined
    # annotation file.  Update the variant row with the count.
    annos = pysam.Tabixfile(args.anno_file)
    naming = guess_contig_naming(annos)
    select_cursor.execute("SELECT chrom, start, end, variant_id FROM variants")
    for row in select_cursor:
        hit_list = []
        for hit in annotations_in_region(row, annos, "tuple", naming=naming):
            try:
                hit_list.append(hit[int(args.col_extract) - 1])
            except IndexError:
                sys.exit("Column " + args.col_extract + " exceeds \
                          the number of columns in your \
                          annotation file. Exiting.")
                          
        hits = ",".join(hit_list)
        if len(hit_list):
            update_qry = "UPDATE variants SET " + args.col_name + " = '" + hits + \
                        "' WHERE variant_id = " + str(row['variant_id'])
        else:
            update_qry = "UPDATE variants SET " + args.col_name + " = NULL" + \
                        " WHERE variant_id = " + str(row['variant_id'])
        update_cursor.execute(update_qry)


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
    select_cursor = conn.cursor()
    update_cursor = conn.cursor()

    if args.col_type == "boolean":
        annotate_variants_bool(args, select_cursor, update_cursor)
    elif args.col_type == "count":
        annotate_variants_count(args, select_cursor, update_cursor)
    elif args.col_type == "list":
        if args.col_extract is None:
            sys.exit("You must specify which column to extract (-e) from the annotation file.")
        else:
            annotate_variants_list(args, select_cursor, update_cursor)
    else:
        sys.exit("Unknown column type requested. Exiting.")


if __name__ == "__main__":
    main()
