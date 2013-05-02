#!/usr/bin/env python

import os
import sys
import sqlite3
import pysam

from gemini.annotations import annotations_in_region, guess_contig_naming

def add_requested_column(col_name, update_cursor):
    """
    Attempt to add a new, user-defined column to the
    variants table.  Exit if the column already exists.
    """
    try:
        alter_qry = "ALTER TABLE variants ADD COLUMN " \
                    + col_name \
                    + " BOOL DEFAULT NULL"
        update_cursor.execute(alter_qry)
    except sqlite3.OperationalError:
        sys.stderr.write("WARNING: column \"("
                         + col_name
                         + ")\" already exists in variants table. Overwriting values.\n")


def _annotate_variants(args, conn, get_val_fn):
    """Generalized annotation of variants with a new column.

    get_val_fn takes a list of annotations in a region and returns
    the value for that region to update the database with.

    Separates selection and identification of values from update,
    to avoid concurrent database access errors from sqlite3, especially on
    NFS systems. The retained to_update list is small, but batching
    could help if memory issues emerge.
    """
    # For each, use Tabix to detect overlaps with the user-defined
    # annotation file.  Update the variant row with T/F if overlaps found.
    anno = pysam.Tabixfile(args.anno_file)
    naming = guess_contig_naming(anno)
    select_cursor = conn.cursor()
    update_cursor = conn.cursor()
    add_requested_column(args.col_name, select_cursor)
    
    last_id = 0
    current_id = 0
    total = 0
    CHUNK_SIZE = 100000
    to_update = []

    select_cursor.execute('''SELECT chrom, start, end, variant_id FROM variants''')
    while True:
        for row in select_cursor.fetchmany(CHUNK_SIZE):
            to_update.append((get_val_fn(annotations_in_region(row,
                                                               anno,
                                                               "tuple",
                                                               naming)),
                              str(row["variant_id"])))
            current_id = row["variant_id"]

        if current_id <= last_id:
            break
        else:
            update_cursor.execute("BEGIN TRANSACTION")
            _update_variants(to_update, args.col_name, update_cursor)
            update_cursor.execute("END TRANSACTION")

            total += len(to_update)
            print "updated", total, "variants"
            last_id = current_id
        to_update = []

def _update_variants(to_update, col_name, cursor):
        update_qry = "UPDATE variants SET " \
                     + col_name \
                     + " = ?" \
                     + " WHERE variant_id = ?"
        cursor.executemany(update_qry, to_update)


def annotate_variants_bool(args, conn):
    """
    Populate a new, user-defined column in the variants
    table with a BOOLEAN indicating whether or not
    overlaps were detected between the variant and the
    annotation file.
    """
    def has_anno_hit(hits):
        for hit in hits:
            return 1
        return 0

    return _annotate_variants(args, conn, has_anno_hit)


def annotate_variants_count(args, conn):
    """
    Populate a new, user-defined column in the variants
    table with a INTEGER indicating the count of overlaps
    between the variant and the
    annotation file.
    """
    def get_hit_count(hits):
        return len(hits)

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
    conn.row_factory = sqlite3.Row  # allow us to refer to columns by name
    conn.isolation_level = None

    if args.col_type == "boolean":
        annotate_variants_bool(args, conn)
    elif args.col_type == "count":
        annotate_variants_count(args, conn)
    elif args.col_type == "list":
        if args.col_extract is None:
            sys.exit("You must specify which column to extract (-e) \
                      from the annotation file.")
        else:
            annotate_variants_list(args, conn)
    else:
        sys.exit("Unknown column type requested. Exiting.")
