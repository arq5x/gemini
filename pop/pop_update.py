#!/usr/bin/env python

import os
import sys
import sqlite3
import pop_constants


def validate_file_for_variants_table(file):
    """
    Ensure that the provided data file is appropriately
    structured for the variant table.  That is, the
    first columns must be the variant_id, and the second
    columns must be the column to be updated or appended.
    """
    fh = open(file, 'r')
    first_line = fh.readline().strip()
    second_line = fh.readline().strip()

    if first_line.startswith("#") and second_line.startswith("#"):
        sys.stderr.write("Data file must have only one header line.\n")
        return False
    elif first_line.startswith("#") and not second_line.startswith("#"):
        fields = first_line[1:].split("\t")
        if len(fields) > 2:
            sys.stderr.write("Variants data file must have two columns: variant_id and user_column.\n")
            return False
        else:
            if fields[0] != pop_constants.VARIANTS_KEY:
                sys.stderr.write("First columns in variants data file must be: %s\n", pop_constants.VARIANTS_KEY)
                return False
            elif len(fields[1]) > 0:
                fh.close()
                return True


def update_table(args, cursor):
    """
    Update a database table based on data provided
    in an external file.
    """
    if args.table == "variants":
        if validate_file_for_variants_table(args.file):
            fh = open(file, 'r')
            first_line = fh.readline().strip()
            fields = first_line[1:].split("\t")
            print fields[1], fields[2]


def alter_table(args, cursor):
    """
    Update a database table based on data provided
    in an external file.
    """
    if args.table == "variants":
        if validate_file_for_variants_table(args.file):
            fh = open(args.file, 'r')
            first_line = fh.readline().strip()
            fields = first_line[1:].split("\t")

            # add the new column to the table
            new_column_name, new_column_type = fields[1].split(":")
            try:
                statement = "alter table variants add column %s %s if not exists default null" % \
                            (new_column_name, new_column_type)
                cursor.execute(statement)
            except sqlite3.Error, e:
                print "An error occurred:", e.args[0]

            try:
                for line in fh:
                    fields = line.strip().split("\t")
                    variant_id = int(fields[0])
                    new_column_value = fields[0]
                    statement = "update variants set %s = %s where variant_id = %s" % \
                                (new_column_name, new_column_value, variant_id)
                    cursor.execute(statement)
            except sqlite3.Error, e:
                print "An error occurred:", e.args[0]
            # load the data into the table



def apply(parser, args):

    if (args.db is None):
        parser.print_help()
        exit()

    if not os.path.exists(args.db):
        sys.stderr.write("Error: cannot find database file.")
        exit(1)
    if not os.path.exists(args.file):
        sys.stderr.write("Error: cannot find update file.")
        exit(1)

    if args.db is not None and args.table is not None:
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        cursor = conn.cursor()
        if args.mode == "update":
            update_table(args, cursor)
        elif args.mode == "append":
            alter_table(args, cursor)


if __name__ == "__main__":
    main()
