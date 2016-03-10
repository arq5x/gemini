#!/usr/bin/env python
import database


def get_table_info(metadata, table_name, out_template):
    """
    Report the column names and types for a given database table
    """
    if not table_name in metadata.tables:
        return
    tbl = metadata.tables[table_name]
    for row in tbl.columns:
        rec = (table_name, row.name, str(row.type))
        print out_template.format(*rec)

def db_info(parser, args):

    conn, metadata = database.get_session_metadata(args.db)
    # column widths for the output
    out_template = "{0:20}{1:30}{2:10}"

        # header
    print out_template.format("table_name", "column_name", "type")
    for table in ['variants', 'variant_impacts', 'samples', 'gene_detailed', 'gene_summary']:
        get_table_info(metadata, table, out_template)
