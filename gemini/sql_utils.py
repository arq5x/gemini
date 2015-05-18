"""
these are utilities to parse and transform SQL statements
"""

import re


def get_select_cols_and_rest(query):
    """
    Separate the a list of selected columns from
    the rest of the query

    Returns:
        1. a list of the selected columns
        2. a string of the rest of the query after the SELECT
    """
    from_loc = query.lower().find("from")

    raw_select_clause = query[0:from_loc].rstrip()
    rest_of_query = query[from_loc:len(query)]

    # remove the SELECT keyword from the query
    select_pattern = re.compile("select", re.IGNORECASE)
    raw_select_clause = select_pattern.sub('', raw_select_clause)

    # now create and iterate through a list of of the SELECT'ed columns
    selected_columns = raw_select_clause.split(',')
    selected_columns = [c.strip() for c in selected_columns]

    return selected_columns, rest_of_query


def ensure_columns(query, cols):
    """
    if a query is missing any of these list of columns, add them
    and return the new query string
    """
    sel_cols, rest = get_select_cols_and_rest(query)
    sel_cols = [x for x in sel_cols]
    if "*" in sel_cols:
        return query
    for c in cols:
        if c not in sel_cols:
            sel_cols += [c]

    sel_string = ", ".join(sel_cols)
    return "select {sel_string} {rest}".format(**locals())
