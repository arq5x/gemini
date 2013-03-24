#!/usr/bin/env python
import sqlite3
import os
import gemini_utils as util


def get_table_info(c, table_name, out_template):
    """
    Report the column names and types for a given database table
    """
    query = "PRAGMA table_info('" + table_name + "')"
    c.execute(query)

    for row in c:
        rec = (table_name, str(row['name']), str(row['type']))
        print out_template.format(*rec)


def db_info(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        # column widths for the output
        out_template = "{0:20}{1:30}{2:10}"

        # header
        print out_template.format("table_name", "column_name", "type")
        for table in ['variants', 'variant_impacts', 'samples']:
            get_table_info(c, table, out_template)
