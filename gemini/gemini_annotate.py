    #!/usr/bin/env python

import os
import sys
import sqlite3
from collections import defaultdict, namedtuple
import atexit
import json
import subprocess
import tempfile

import numpy as np
from scipy.stats import mode
import pysam
import cyvcf as vcf

from gemini.annotations import annotations_in_region, annotations_in_vcf, guess_contig_naming
from database import database_transaction

def add_requested_columns(args, update_cursor, col_names, col_types=None):
    """
    Attempt to add new, user-defined columns to the
    variants table.  Warn if the column already exists.
    """

    if args.anno_type in ["count", "boolean"]:

        col_name = col_names[0]
        col_type = "integer"

        try:
            alter_qry = "ALTER TABLE variants ADD COLUMN " \
                        + col_name \
                        + " " \
                        + col_type \
                        + " " \
                        + "DEFAULT NULL"
            update_cursor.execute(alter_qry)
        except sqlite3.OperationalError:
            sys.stderr.write("WARNING: Column \"("
                             + col_name
                             + ")\" already exists in variants table. Overwriting values.\n")
            # reset values so that records don't retain old annotations.
            update_cursor.execute("UPDATE variants SET " + col_name + " = NULL WHERE 1")

    elif args.anno_type == "extract":

        for col_name, col_type in zip(col_names, col_types):

            try:
                alter_qry = "ALTER TABLE variants ADD COLUMN " \
                            + col_name \
                            + " " \
                            + col_type \
                            + " " \
                            + "DEFAULT NULL"
                update_cursor.execute(alter_qry)
            except sqlite3.OperationalError:
                sys.stderr.write("WARNING: Column \"("
                                 + col_name
                                 + ")\" already exists in variants table. Overwriting values.\n")
    else:
        sys.exit("Unknown annotation type: %s\n" % args.anno_type)


def _annotate_variants(args, conn, get_val_fn, col_names=None, col_types=None, col_ops=None):
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
    add_requested_columns(args, select_cursor, col_names, col_types)

    last_id = 0
    current_id = 0
    total = 0
    CHUNK_SIZE = 100000
    to_update = []

    select_cursor.execute('''SELECT chrom, start, end, ref, alt, variant_id FROM variants''')
    while True:
        for row in select_cursor.fetchmany(CHUNK_SIZE):

            # update_data starts out as a list of the values that should
            # be used to populate the new columns for the current row.
            # Prefer no pysam parsing over tuple parsing to work around bug in pysam 0.8.0
            # https://github.com/pysam-developers/pysam/pull/44
            if args.anno_file.endswith(('.vcf', '.vcf.gz')):
                update_data = get_val_fn(annotations_in_vcf(row, anno, None, naming, args.region_only, True))
            else:
                update_data = get_val_fn(annotations_in_region(row, anno, None, naming))
            #update_data = get_val_fn(annotations_in_region(row, anno, "tuple", naming))
            # were there any hits for this row?
            if len(update_data) > 0:
                # we add the primary key to update_data for the
                # where clause in the SQL UPDATE statement.
                update_data.append(str(row["variant_id"]))
                to_update.append(tuple(update_data))

            current_id = row["variant_id"]

        if current_id <= last_id:
            break
        else:
            update_cursor.execute("BEGIN TRANSACTION")
            _update_variants(to_update, col_names, update_cursor)
            update_cursor.execute("END TRANSACTION")

            total += len(to_update)
            print "updated", total, "variants"
            last_id = current_id
        to_update = []

def _update_variants(to_update, col_names, cursor):
    update_qry = "UPDATE variants SET "
    update_cols = ",".join(col_name + " = ?" for col_name in col_names)
    update_qry += update_cols
    update_qry += " WHERE variant_id = ?"
    cursor.executemany(update_qry, to_update)

def annotate_variants_bool(args, conn, col_names):
    """
    Populate a new, user-defined column in the variants
    table with a BOOLEAN indicating whether or not
    overlaps were detected between the variant and the
    annotation file.
    """
    def has_hit(hits):
        for hit in hits:
            return [1]
        return [0]

    return _annotate_variants(args, conn, has_hit, col_names)


def annotate_variants_count(args, conn, col_names):
    """
    Populate a new, user-defined column in the variants
    table with a INTEGER indicating the count of overlaps
    between the variant and the
    annotation file.
    """
    def get_hit_count(hits):
        return [len(list(hits))]

    return _annotate_variants(args, conn, get_hit_count, col_names)


def _map_list_types(hit_list, col_type):
    # TODO: handle missing because of VCF.
    try:
        if col_type == "int":
            return [int(h) for h in hit_list if not h in (None, 'nan')]
        elif col_type == "float":
            return [float(h) for h in hit_list if not h in (None, 'nan')]
    except ValueError:
        sys.exit('Non-numeric value found in annotation file: %s\n' % (','.join(hit_list)))


def gemops_mean(li, col_type):
    return np.average(_map_list_types(li, col_type))

def gemops_sum(li, col_type):
    return np.sum(_map_list_types(li, col_type))

def gemops_list(li, col_type):
    return ",".join(li)

def gemops_uniq_list(li, col_type):
    return ",".join(set(li))

def gemops_median(li, col_type):
    return np.median(_map_list_types(li, col_type))

def gemops_min(li, col_type):
    return np.min(_map_list_types(li, col_type))

def gemops_max(li, col_type):
    return np.max(_map_list_types(li, col_type))

def gemops_mode(li, col_type):
    return mode(_map_list_types(li, col_type))[0][0]

def gemops_first(li, col_type):
    return li[0]

def gemops_last(li, col_type):
    return li[-1]

# lookup from the name to the func above.
op_funcs = dict((k[7:], v) for k, v in locals().items() if k.startswith('gemops_'))

def fix_val(val, type):
    if not type in ("int", "float"): return val
    if isinstance(val, (int, float)): return val

    if type == "int": fn = int
    else: fn = float
    if not val:
        return None
    try:
        return fn(val)
    except ValueError:
        sys.exit('Non %s value found in annotation file: %s\n' % (type, val))

def get_hit_list(hits, col_idxs, args):
    hits = list(hits)
    if len(hits) == 0:
        return []

    hit_list = defaultdict(list)
    for hit in hits:
        if isinstance(hit, basestring):
            hit = hit.split("\t")
        if args.anno_file.endswith(('.vcf', '.vcf.gz')):
            # only makes sens to extract when there is an equal sign
            info = dict((x[0], x[1]) for x in (p.split('=') for p in hit[7].split(';') if '=' in p))
            for idx, col_idx in enumerate(col_idxs):
                if not col_idx in info:
                    hit_list[idx].append('nan')
                    sys.stderr.write("WARNING: %s is missing from INFO field in %s for at "
                        "least one record.\n" % (col_idx, args.anno_file))
                else:
                    hit_list[idx].append(info[col_idx])
                # just append None since in a VCFthey are likely # to be missing ?


        else:
            try:
                for idx, col_idx in enumerate(col_idxs):
                    hit_list[idx].append(hit[int(col_idx) - 1])
            except IndexError:
                sys.exit("EXITING: Column " + args.col_extracts + " exceeds "
                          "the number of columns in your "
                          "annotation file.\n")
    return hit_list

def annotate_variants_extract(args, conn, col_names, col_types, col_ops, col_idxs):
    """
    Populate a new, user-defined column in the variants
    table based on the value(s) from a specific column.
    in the annotation file.
    """
    def summarize_hits(hits):
        hit_list = get_hit_list(hits, col_idxs, args)
        if hit_list == []: return []
        vals = []
        for idx, op in enumerate(col_ops):
            # more than one overlap, must summarize
            try:
                val = op_funcs[op](hit_list[idx], col_types[idx])
            except ValueError:
                val = None
            vals.append(fix_val(val, col_types[idx]))
        return vals

    return _annotate_variants(args, conn, summarize_hits,
                              col_names, col_types, col_ops)

def annotate(parser, args):
    def _validate_args(args):
        if (args.col_operations or args.col_types or args.col_extracts):
            sys.exit('EXITING: You may only specify a column name (-c) when '
                     'using \"-a boolean\" or \"-a count\".\n')

        col_names = args.col_names.split(',')
        if len(col_names) > 1:
            sys.exit('EXITING: You may only specify a single column name (-c) '
                     'when using \"-a boolean\" or \"-a count\".\n')

        if not args.anno_file.endswith(('.vcf', '.vcf.gz')) and args.region_only and parser is not None:
            sys.exit('EXITING: You may only specify --region-only when annotation is a VCF.')

        return col_names

    def _validate_extract_args(args):
        if args.anno_file.endswith(('.vcf', '.vcf.gz')):
            if not args.col_names:
                args.col_names = args.col_extracts
            elif not args.col_extracts:
                args.col_extracts = args.col_names
        elif args.region_only and parser is not None:
            sys.exit('EXITING: You may only specify --region-only when annotation is a VCF.1')

        if not args.col_types:
            sys.exit('EXITING: need to give column types ("-t")\n')
        col_ops = args.col_operations.split(',')
        col_idxs = args.col_extracts.split(',')

        col_names = args.col_names.split(',')
        col_types = args.col_types.split(',')

        supported_types = ['text', 'float', 'integer']
        for col_type in col_types:
            if col_type not in supported_types:
                sys.exit('EXITING: Column type [%s] not supported.\n' %
                         (col_type))

        supported_ops = op_funcs.keys()

        for col_op in col_ops:
            if col_op not in supported_ops:
                sys.exit('EXITING: Column operation [%s] not supported.\n' %
                         (col_op))

        if not (len(col_ops) == len(col_names) ==
                len(col_types) == len(col_idxs)):
            sys.exit('EXITING: The number of column names, numbers, types, and '
                     'operations must match: [%s], [%s], [%s], [%s]\n' %
                     (args.col_names, args.col_extracts, args.col_types, args.col_operations))

        return col_names, col_types, col_ops, col_idxs

    if (args.db is None):
        parser.print_help()
        exit(1)
    if not os.path.exists(args.db):
        sys.stderr.write("Error: cannot find database file.")
        exit(1)
    if not os.path.exists(args.anno_file):
        sys.stderr.write("Error: cannot find annotation file.")
        exit(1)

    conn = sqlite3.connect(args.db)
    conn.row_factory = sqlite3.Row  # allow us to refer to columns by name
    conn.isolation_level = None

    if args.anno_type == "boolean":
        col_names = _validate_args(args)
        annotate_variants_bool(args, conn, col_names)
    elif args.anno_type == "count":
        col_names = _validate_args(args)
        annotate_variants_count(args, conn, col_names)
    elif args.anno_type == "extract":
        if args.col_extracts is None and not args.anno_file.endswith('.vcf.gz'):
            sys.exit("You must specify which column to "
                     "extract from your annotation file.")
        else:
            col_names, col_types, col_ops, col_idxs = _validate_extract_args(args)
            annotate_variants_extract(args, conn, col_names, col_types, col_ops, col_idxs)
    else:
        sys.exit("Unknown column type requested. Exiting.")

    conn.close()

    # index on the newly created columns
    for col_name in col_names:
        with database_transaction(args.db) as c:
            c.execute('''drop index if exists %s''' % (col_name + "idx"))
            c.execute('''create index %s on variants(%s)''' % (col_name + "idx", col_name))


# ## Automate addition of extra fields to database

def add_extras(gemini_db, chunk_dbs, region_only):
    """Annotate gemini database with extra columns from processed chunks, if available.
    """
    for chunk in chunk_dbs:
        extra_file = get_extra_vcf(chunk)
        if extra_file is False:
            # there was not extra annotation so we just continue
            continue
        # these the the field names that we'll pull from the info field.
        fields = [x.strip() for x in open(extra_file[:-3] + ".fields")]

        ops = ["first" for t in fields]
        Args = namedtuple("Args", "db,anno_file,anno_type,col_operations,col_names,col_types,col_extracts,region_only")

        # TODO: hard-coded "text" into the type...
        args = Args(gemini_db, extra_file, "extract", ",".join(ops),
                    ",".join(fields), ",".join(["text"] * len(fields)),
                    ",".join(fields),
                    region_only)
        annotate(None, args)
        os.unlink(extra_file[:-3] + ".fields")


def rm(path):
    try:
        os.unlink(path)
    except:
        pass


def get_extra_vcf(gemini_db, tmpl=None):
    """Retrieve extra file associated with a gemini database.
    Most commonly, this will be with VEP annotations added.
    Returns false if there are no vcfs associated with the database.
    """
    base = os.path.basename(gemini_db)
    path = os.path.join(tempfile.gettempdir(), "extra.%s.vcf" % base)
    mode = "r" if tmpl is None else "w"
    if mode == "r":
        if not os.path.exists(path):
            return False
        if not path.endswith(".gz"):
            subprocess.check_call(["bgzip", "-f", path])
            bgzip_out = path + ".gz"
            subprocess.check_call(["tabix", "-p", "vcf", "-f", bgzip_out])
            return bgzip_out
        return path

    fh = open(path, "w")
    if mode == "w":
        atexit.register(rm, fh.name)
        atexit.register(rm, fh.name + ".gz")
        atexit.register(rm, fh.name + ".gz.tbi")

        return vcf.Writer(fh, tmpl)
    return vcf.Reader(fh)
