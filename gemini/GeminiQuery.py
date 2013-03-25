#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import string
import itertools
from pyparsing import ParseResults
from collections import defaultdict
#from collections import OrderedDict

# gemini imports
import gemini_utils as util
from gemini_constants import *
from gemini_utils import OrderedSet, OrderedDict
import compression


class GeminiRow(object):

    def __init__(self, row):
        #super(OrderedDict, self).__init__()
        self.row = row

    def __getitem__(self, val):
        return self.row[val]

    def __repr__(self):
        return '\t'.join([str(self.row[c]) for c in self.row])


class GeminiQuery(object):

    def __init__(self, db):
        self.db = db
        self.query_executed = False
        self.connected_to_db = False

    def run(self, query, gt_filter=None, use_header=False):
        """
        Execute the query against the Gemini database
        """
        self.query = query
        self.gt_filter = gt_filter
        self.use_header = use_header
        self.header_processed = False

        if not self.connected_to_db:
            self._connect_to_database()
            self.sample_to_idx = util.map_samples_to_indicies(self.c)
            self.connected_to_db = True

        self.query_pieces = self.query.split()
        if not any(s.startswith("gt") for s in self.query_pieces) and \
                not any("gt" in s for s in self.query_pieces):
            self.query_type = "no-genotypes"
            self._apply_basic_query()
        else:
            if self.gt_filter is None:
                self.query_type = "select-genotypes"
            # using a gt-filter
            else:
                self.gt_filter = self._correct_genotype_filter()
                self.query_type = "filter-genotypes"
            self._apply_query_w_genotypes()

        self.query_executed = True

    def __iter__(self):
        return self
    
    def next(self):
        """
        Return a GeminiRow object for the next query result
        """

        # handle a header if requested
        if self.use_header and not self.header_processed:
            self.header_processed = True

            if self.query_type == "no-genotypes":
                h = [col for col in self.all_query_cols]
                return GeminiRow(OrderedDict(itertools.izip(h, h)))
            else:
                h = [col for col in self.all_query_cols] + \
                    [col for col in OrderedSet(self.all_cols_orig) - OrderedSet(self.select_cols)]
                return GeminiRow(OrderedDict(itertools.izip(h, h)))

        try:
            row = self.c.next()

            if self.query_type == "no-genotypes":
                return GeminiRow(OrderedDict(row))
            # handle queries that work with genotypes
            else:
                gts = compression.unpack_genotype_blob(row['gts'])
                gt_types = compression.unpack_genotype_blob(row['gt_types'])
                gt_phases = compression.unpack_genotype_blob(row['gt_phases'])
                gt_depths = compression.unpack_genotype_blob(row['gt_depths'])
                
                # skip the row if it fails the requested genotype filter
                if self.query_type == "filter-genotypes" and not eval(self.gt_filter):
                    return self.next()

                fields = OrderedDict()
                for idx, col in enumerate(self.report_cols):
                    if col == "*":
                        continue
                    if not col.startswith("gt") and not col.startswith("GT"):
                        fields[col] = row[col]
                    else:
                        # reuse the original column anme user requested
                        # e.g. replace gts[1085] with gts.NA20814
                        orig_col = self.all_cols_orig[idx]
                        fields[orig_col] = eval(col.strip())
                return GeminiRow(fields)
        except:
            raise StopIteration

    def _connect_to_database(self):
        """
        Establish a connection to the requested Gemini database
        """
        # open up a new database
        if os.path.exists(self.db):
            self.conn = sqlite3.connect(self.db)
            self.conn.isolation_level = None
            # allow us to refer to columns by name
            self.conn.row_factory = sqlite3.Row 
            self.c = self.conn.cursor()

    def _apply_basic_query(self):
        """
        Execute a vanilla query. That is, not gt* columns
        are in either the select or where clauses.
        """
        self.c.execute(self.query)
        self.all_query_cols = [str(tuple[0]) for tuple in self.c.description
                    if not tuple[0].startswith("gt")]

    def _apply_query_w_genotypes(self):
        """
        Execute a query that contains gt* columns in only in the SELECT.
        """
    
        (self.select_cols, self.all_cols_new, 
         self.all_cols_orig, self.gt_col_map) = self._split_select()
    
        self.query = self._add_gt_cols_to_query()
        self.c.execute(self.query)
    
        # what are the columns that were actually selected by the user.
        self.all_query_cols = [str(tuple[0]) for tuple in self.c.description
                          if not tuple[0].startswith("gt")]
    
        if "*" in self.select_cols:
            self.select_cols.remove("*")
            all_cols_orig.remove("*")
            self.all_cols_new.remove("*")
            self.select_cols += self.all_query_cols

        self.report_cols = self.all_query_cols + \
                            list(OrderedSet(self.all_cols_new) - OrderedSet(self.select_cols))

    def _correct_genotype_col(self, raw_col):
        """
        Convert a _named_ genotype index to a _numerical_
        genotype index so that the appropriate value can be
        extracted for the sample from the genotype numpy arrays.
        
        These lookups will be eval()'ed on the resuting rows to
        extract the appropriate information.
        
        For example, convert gt_types.1478PC0011 to gt_types[11]
        """
        if raw_col == "*":
            return raw_col.lower()
        (column, sample) = raw_col.split('.')
        corrected = column.lower() + "[" + str(self.sample_to_idx[sample]).lower() + "]"
        return corrected

    def _correct_genotype_filter(self):
        """
        This converts a "raw" genotype filter supplied by the user
        to a filter than can be eval()'ed.  Specifically, we must
        convery a _named_ genotype index to a _numerical_
        genotype index so that the appropriate value can be
        extracted for the sample from the genotype numpy arrays.

        For example, converts:
        --gt-filter "(gt_types.1478PC0011 == 1)"
        to
        (gt_types[11] == 1)
        """
        corrected_gt_filter = []
        tokens = re.split(r'[\s+]+', self.gt_filter)
        for token in tokens:
            if token.find("gt") >= 0 or token.find("GT") >= 0:
                corrected = self._correct_genotype_col(token)
                corrected_gt_filter.append(corrected)
            else:
                corrected_gt_filter.append(token)
        return " ".join(corrected_gt_filter)

    def _add_gt_cols_to_query(self):
        """
        We have to modify the raw query to select the genotype
        columns in order to support the genotype filters.  That is,
        if the user wants to limit the rows returned based upon, for example,
        "gts.joe == 1", then we need to select the full gts BLOB column in
        order to enforce that limit.  The user wouldn't have selected gts as a
        columns, so therefore, we have to modify the select statement to add
        it.
        
        In essence, when a gneotype filter has been requested, we always add
        the gts, gt_types and gt_phases columns.
        """
        from_loc = self.query.lower().find("from")
        if from_loc > 1:
            raw_select_clause = self.query[0:from_loc].rstrip()
            rest_of_query = self.query[from_loc:len(self.query)]
        
            # remove any GT columns
            select_clause_list = []
            for token in raw_select_clause.split():
                if not token.startswith("gt") and not token.startswith("GT"):
                    select_clause_list.append(token)
        
            # add the genotype columns to the query
            if select_clause_list[len(select_clause_list) - 1].endswith(",") or \
                (len(select_clause_list) == 1 and
                 select_clause_list[0].strip().lower() == "select"):
                select_clause = " ".join(select_clause_list) + \
                    " gts, gt_types, gt_phases, gt_depths "
            else:
                select_clause = " ".join(select_clause_list) + \
                    ", gts, gt_types, gt_phases, gt_depths "
            self.query = select_clause + rest_of_query
            # extract the original select columns
            return self.query
        
        else:
            sys.exit("Malformed query.")
    
    def _split_select(self):
        """
        Build a list of _all_ columns in the SELECT statement
        and segregated the non-genotype specific SELECT columns.
        
        This is used to control how to report the results, as the
        genotype-specific columns need to be eval()'ed whereas others
        do not.
        
        For example: "SELECT chrom, start, end, gt_types.1478PC0011"
        will populate the lists as follows:
        
        select_columns = ['chrom', 'start', 'end']
        all_columns = ['chrom', 'start', 'end', 'gt_types[11]']
        """
        select_columns = []
        all_columns_new = []
        all_columns_orig = []
        gt_col_map = {}
        
        # iterate through all of the select columns andclear
        # distinguish the genotype-specific columns from the base columns
        from_loc = self.query.lower().find("from")
        if from_loc < 1:
            sys.exit("Malformed query.")
        
        raw_select_clause = self.query[0:from_loc].rstrip()
        rest_of_query = self.query[from_loc:len(query)]
        
        for token in raw_select_clause.replace(',', '').split():
            if token == "SELECT" or token == "select":
                continue
            if not token.startswith("GT") and not token.startswith("gt"):
                select_columns.append(token)
                all_columns_new.append(token)
                all_columns_orig.append(token)
            else:
                new_col = self._correct_genotype_col(token)
                all_columns_new.append(new_col)
                all_columns_orig.append(token)
                gt_col_map[token] = new_col
        
        return select_columns, all_columns_new, all_columns_orig, gt_col_map


if __name__ == "__main__":

    db = sys.argv[1]

    g_query = GeminiQuery(db)

    print "test a basic query with no genotypes"
    query  = "select chrom, start, end from variants limit 5"
    g_query.run(query)
    for row in g_query:
        print row


    print "test a basic query with no genotypes using a header"
    query  = "select chrom, start, end from variants limit 5"
    g_query.run(query, use_header=True)
    for row in g_query:
        print row


    print "test query that selects a sample genotype"
    query  = "select chrom, start, end, gts.NA20814 from variants limit 5"
    g_query.run(query)
    for row in g_query:
        print row


    print "test query that selects a sample genotype and uses a header"
    query  = "select chrom, start, end, gts.NA20814 from variants limit 5"
    g_query.run(query, use_header=True)
    for row in g_query:
        print row


    print "test query that selects and _filters_ on a sample genotype"
    query  = "select chrom, start, end, gts.NA20814 from variants limit 50"
    db_filter = "gt_types.NA20814 == HET"
    g_query.run(query, db_filter)
    for row in g_query:
        print row

    print "test query that selects and _filters_ on a sample genotype and uses a filter"
    query  = "select chrom, start, end, gts.NA20814 from variants limit 50"
    db_filter = "gt_types.NA20814 == HET"
    g_query.run(query, db_filter, use_header=True)
    for row in g_query:
        print row

