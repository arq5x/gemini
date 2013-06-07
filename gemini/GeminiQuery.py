#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import itertools
import collections

# gemini imports
import gemini_utils as util
from gemini_constants import *
from gemini_utils import OrderedSet, OrderedDict
import compression


class GeminiRow(object):

    def __init__(self, row, gts=None, gt_types=None,
                 gt_phases=None, gt_depths=None):
        self.row = row
        self.gts = gts
        self.gt_types = gt_types
        self.gt_phases = gt_phases
        self.gt_depths = gt_depths
        self.gt_cols = ['gts', 'gt_types', 'gt_phases', 'gt_depths']

    def __getitem__(self, val):
        if val not in self.gt_cols:
            return self.row[val]
        else:
            if val == 'gts':
                return self.gts
            elif val == 'gt_types':
                return self.gt_types
            elif val == 'gt_phases':
                return self.gt_phases
            elif val == 'gt_depths':
                return self.gt_depths

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.row.keys()
        except:
            raise StopIteration

    def __repr__(self):
        return '\t'.join([str(self.row[c]) for c in self.row])


class GeminiQuery(object):

    """
    An interface to submit queries to an existing Gemini database
    and iterate over the results of the query.

    We create a GeminiQuery object by specifying database to which to
    connect::
        from gemini import GeminiQuery
        gq = GeminiQuery("my.db")

    We can then issue a query against the database and iterate through
    the results by using the ``run()`` method::

        gq.run("select chrom, start, end from variants")
        for row in gq:
            print row

    Instead of printing the entire row, one access print specific columns::

        gq.run("select chrom, start, end from variants")
        for row in gq:
            print row['chrom']

    Also, all of the underlying numpy genotype arrays are
    always available::

        gq.run("select chrom, start, end from variants")
        for row in gq:
            gts = row.gts
            print row['chrom'], gts
            # yields "chr1" ['A/G' 'G/G' ... 'A/G']

    The ``run()`` methods also accepts genotype filter::

        query = "select chrom, start, end" from variants"
        gt_filter = "gt_types.NA20814 == HET"
        gq.run(query)
        for row in gq:
            print row

    Lastly, one can use the ``sample_to_idx`` and ``idx_to_sample``
    dictionaries to gain access to sample-level genotype information
    either by sample name or by sample index::

        # grab dict mapping sample to genotype array indices
        smp2idx = gq.sample_to_idx

        query  = "select chrom, start, end from variants"
        gt_filter  = "gt_types.NA20814 == HET"
        gq.run(query, gt_filter)

        # print a header listing the selected columns
        print gq.header
        for row in gq:
            # access a NUMPY array of the sample genotypes.
            gts = row['gts']
            # use the smp2idx dict to access sample genotypes
            idx = smp2idx['NA20814']
            print row, gts[idx]
    """

    def __init__(self, db):
        self.db = db
        self.query_executed = False
        self.for_browser = False

        self._connect_to_database()
        # map sample names to indices. e.g. self.sample_to_idx[NA20814] -> 323
        self.sample_to_idx = util.map_samples_to_indicies(self.c)
        # and vice versa. e.g., self.idx_to_sample[323] ->  NA20814
        self.idx_to_sample = util.map_indicies_to_samples(self.c)

    def _set_gemini_browser(self, for_browser):
        self.for_browser = for_browser

    def run(self, query, gt_filter=None, show_variant_samples=False):
        """
        Execute a query against a Gemini database. The user may
        specify:

            1. (reqd.) an SQL `query`.
            2. (opt.) a genotype filter.
        """
        self.query = query
        self.gt_filter = gt_filter
        self.show_variant_samples = show_variant_samples

        self.query_pieces = self.query.split()
        if not any(s.startswith("gt") for s in self.query_pieces) and \
                not any("gt" in s for s in self.query_pieces):
            if self.gt_filter is None:
                self.query_type = "no-genotypes"
            else:
                self.gt_filter = self._correct_genotype_filter()
                self.query_type = "filter-genotypes"
        else:
            if self.gt_filter is None:
                self.query_type = "select-genotypes"
            else:
                self.gt_filter = self._correct_genotype_filter()
                self.query_type = "filter-genotypes"

        self._apply_query()
        self.query_executed = True


    def __iter__(self):
        return self

    @property
    def header(self):
        """
        Return a header describing the columns that
        were selected in the query issued to a GeminiQuery object.
        """

        if self.query_type == "no-genotypes":
            h = [col for col in self.all_query_cols]
        else:
            h = [col for col in self.all_query_cols] + \
                [col for col in OrderedSet(self.all_columns_orig)
                 - OrderedSet(self.select_columns)]
        if self.show_variant_samples:
            h += ["samples_with_variant", "HET_samples", "HOM_ALT_samples"]
        return GeminiRow(OrderedDict(itertools.izip(h, h)))

    @property
    def sample2index(self):
        """
        Return a dictionary mapping sample names to
        genotype array offsets::

            gq = GeminiQuery("my.db")
            s2i = gq.sample2index

            print s2i['NA20814']
            # yields 1088
        """
        return self.sample_to_idx

    @property
    def index2sample(self):
        """
        Return a dictionary mapping sample names to
        genotype array offsets::

            gq = GeminiQuery("my.db")
            i2s = gq.index2sample

            print i2s[1088]
            # yields "NA20814"
        """
        return self.idx_to_sample

    def next(self):
        """
        Return the GeminiRow object for the next query result.
        """
        # we use a while loop since we may skip records based upon
        # genotype filters.  if we need to skip a record, we just
        # throw a continue and keep trying. the alternative is to just
        # recursively call self.next() if we need to skip, but this
        # can quickly exceed the stack.
        while (1):
            try:
                row = self.c.next()

                if self._query_needs_genotype_info():
                    gts = compression.unpack_genotype_blob(row['gts'])
                    gt_types = \
                        compression.unpack_genotype_blob(row['gt_types'])
                    gt_phases = \
                        compression.unpack_genotype_blob(row['gt_phases'])
                    gt_depths = \
                        compression.unpack_genotype_blob(row['gt_depths'])

                    # skip the record if it does not meet the user's genotype filter
                    if self.gt_filter and not eval(self.gt_filter):
                        continue

                fields = OrderedDict()

                for idx, col in enumerate(self.report_cols):
                    if col == "*":
                        continue
                    if not col.startswith("gt") and not col.startswith("GT"):
                        fields[col] = row[col]
                    else:
                        # reuse the original column anme user requested
                        # e.g. replace gts[1085] with gts.NA20814
                        if '[' in col:
                            orig_col = self.gt_idx_to_name_map[col]
                            fields[orig_col] = eval(col.strip())
                        else:
                            # asked for "gts" or "gt_types", e.g.
                            if col == "gts":
                                fields[col] = ','.join(gts)
                            elif col == "gt_types":
                                fields[col] = \
                                    ','.join(str(t) for t in gt_types)
                            elif col == "gt_phases":
                                fields[col] = \
                                    ','.join(str(p) for p in gt_phases)
                            elif col == "gt_depths":
                                fields[col] = \
                                    ','.join(str(d) for d in gt_depths)

                if self.show_variant_samples:
                    gt_types = compression.unpack_genotype_blob(row['gt_types'])
                    variant_samples = [x for x, y in enumerate(gt_types) if y == HET or
                                       y == HOM_ALT]
                    variant_names = [self.idx_to_sample[x] for x in variant_samples]
                    fields["variant_samples"] = ",".join(variant_names)
                    het_samples = [x for x, y in enumerate(gt_types) if y == HET]
                    het_names = [self.idx_to_sample[x] for x in het_samples]
                    fields["HET_samples"] = ",".join(het_names)
                    hom_alt_samples = [x for x, y in enumerate(gt_types) if y == HOM_ALT]
                    hom_alt_names = [self.idx_to_sample[x] for x in hom_alt_samples]
                    fields["HOM_ALT_samples"] = ",".join(hom_alt_names)

                if self._query_needs_genotype_info():
                    if not self.for_browser:
                        return GeminiRow(fields,
                                         gts, gt_types, gt_phases, gt_depths)
                    else:
                        return fields
                else:
                    if not self.for_browser:
                        return GeminiRow(fields)
                    else:
                        return fields
            except:
                raise StopIteration

    def _connect_to_database(self):
        """
        Establish a connection to the requested Gemini database.
        """
        # open up a new database
        if os.path.exists(self.db):
            self.conn = sqlite3.connect(self.db)
            self.conn.isolation_level = None
            # allow us to refer to columns by name
            self.conn.row_factory = sqlite3.Row
            self.c = self.conn.cursor()

    def _apply_query(self):
        """
        Execute a query. Intercept gt* columns and
        replace sample names with indices where necessary.
        """
        if self._query_needs_genotype_info():
            # break up the select statement into individual
            # pieces and replace genotype columns using sample
            # names with sample indices
            self._split_select()

            # we only need genotype information if the user is
            # querying the variants table
            self.query = self._add_gt_cols_to_query()

            self.c.execute(self.query)

            self.all_query_cols = [str(tuple[0]) for tuple in self.c.description
                                   if not tuple[0].startswith("gt")]

            if "*" in self.select_columns:
                self.select_columns.remove("*")
                self.all_columns_orig.remove("*")
                self.all_columns_new.remove("*")
                self.select_columns += self.all_query_cols

            self.report_cols = self.all_query_cols + \
                list(OrderedSet(self.all_columns_new) - OrderedSet(self.select_columns))
        # the query does not involve the variants table
        # and as such, we don't need to do anything fancy.
        else:
            self.c.execute(self.query)
            self.all_query_cols = [str(tuple[0]) for tuple in self.c.description
                                   if not tuple[0].startswith("gt")]
            self.report_cols = self.all_query_cols

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
        elif '.' in raw_col:
            # e.g., "gts.NA12878"
            (column, sample) = raw_col.split('.')
            corrected = column.lower() + "[" + str(self.sample_to_idx[sample]).lower() + "]"
        else:
            # e.g. "gts" - do nothing
            corrected = raw_col
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
        self.select_columns = []
        self.all_columns_new = []
        self.all_columns_orig = []
        self.gt_name_to_idx_map = {}
        self.gt_idx_to_name_map = {}

        # iterate through all of the select columns andclear
        # distinguish the genotype-specific columns from the base columns
        from_loc = self.query.lower().find("from")
        if from_loc < 1:
            sys.exit("Malformed query.")

        raw_select_clause = self.query[0:from_loc].rstrip()
        rest_of_query = self.query[from_loc:len(self.query)]

        for token in raw_select_clause.replace(',', '').split():
            if token == "SELECT" or token == "select":
                continue
            if not token.startswith("GT") and not token.startswith("gt"):
                self.select_columns.append(token)
                self.all_columns_new.append(token)
                self.all_columns_orig.append(token)
            else:
                new_col = self._correct_genotype_col(token)
                self.all_columns_new.append(new_col)
                self.all_columns_orig.append(token)
                self.gt_name_to_idx_map[token] = new_col
                self.gt_idx_to_name_map[new_col] = token

    def _tokenize_query(self):
        tokens = list(flatten([x.split(",") for x in self.query.split(" ")]))
        return tokens

    def _query_needs_genotype_info(self):
        tokens = self._tokenize_query()
        requested_genotype = "variants" in tokens and any([x.startswith("gt") for x in tokens])
        return requested_genotype or self.show_variant_samples

def flatten(l):
    """
    flatten an irregular list of lists
    example: flatten([[[1, 2, 3], [4, 5]], 6]) -> [1, 2, 3, 4, 5, 6]
    lifted from: http://stackoverflow.com/questions/2158395/

    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                   basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

if __name__ == "__main__":

    db = sys.argv[1]

    gq = GeminiQuery(db)

    print "test a basic query with no genotypes"
    query = "select chrom, start, end from variants limit 5"
    gq.run(query)
    for row in gq:
        print row

    print "test a basic query with no genotypes using a header"
    query = "select chrom, start, end from variants limit 5"
    gq.run(query)
    print gq.header
    for row in gq:
        print row

    print "test query that selects a sample genotype"
    query = "select chrom, start, end, gts.NA20814 from variants limit 5"
    gq.run(query)
    for row in gq:
        print row

    print "test query that selects a sample genotype and uses a header"
    query = "select chrom, start, end, gts.NA20814 from variants limit 5"
    gq.run(query)
    print gq.header
    for row in gq:
        print row

    print "test query that selects and _filters_ on a sample genotype"
    query = "select chrom, start, end, gts.NA20814 from variants limit 50"
    db_filter = "gt_types.NA20814 == HET"
    gq.run(query, db_filter)
    for row in gq:
        print row

    print "test query that selects and _filters_ on a sample genotype and uses a filter"
    query = "select chrom, start, end, gts.NA20814 from variants limit 50"
    db_filter = "gt_types.NA20814 == HET"
    gq.run(query, db_filter)
    print gq.header
    for row in gq:
        print row

    print "test query that selects and _filters_ on a sample genotype and uses a filter and a header"
    query = "select chrom, start, end, gts.NA20814 from variants limit 50"
    db_filter = "gt_types.NA20814 == HET"
    gq.run(query, db_filter)
    print gq.header
    for row in gq:
        print row

    print "demonstrate accessing individual columns"
    query = "select chrom, start, end, gts.NA20814 from variants limit 50"
    db_filter = "gt_types.NA20814 == HET"
    gq.run(query, db_filter)
    for row in gq:
        print row['chrom'], row['start'], row['end'], row['gts.NA20814']
