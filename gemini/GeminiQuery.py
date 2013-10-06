#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import itertools
import collections
import json
import abc
import re

# gemini imports
import gemini_utils as util
from gemini_constants import *
from gemini_utils import OrderedSet, OrderedDict, itersubclasses
import compression
from sql_utils import ensure_columns, get_select_cols_and_rest


class RowFormat:
    """A row formatter to output rows in a custom format.  To provide
    a new output format 'foo', implement the class methods and set the
    name field to foo.  This will automatically add support for 'foo' to
    anything accepting the --format option via --format foo.
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def name(self):
        return

    @abc.abstractmethod
    def format(self, row):
        """ return a string representation of a GeminiRow object
        """
        return '\t'.join([str(row[c]) for c in row])

    @abc.abstractmethod
    def format_query(self, query):
        """ augment the query with columns necessary for the format or else just
        return the untouched query
        """
        return query

    @abc.abstractmethod
    def predicate(self, row):
        """ the row must pass this additional predicate to be output. Just
        return True if there is no additional predicate"""
        return True


class DefaultRowFormat(RowFormat):

    name = "default"

    @classmethod
    def format(self, row):
        return '\t'.join([str(row[c]) for c in row])

    @classmethod
    def format_query(self, query):
        return query

    @classmethod
    def predicate(self, row):
        return True


class TPEDRowFormat(RowFormat):

    name = "tped"
    NULL_GENOTYPES = ["."]

    @classmethod
    def format(self, row):
        VALID_CHROMOSOMES = map(str, range(1, 23)) + ["X", "Y", "XY", "MT"]
        chrom = row['chrom'].split("chr")[1]
        chrom = chrom if chrom in VALID_CHROMOSOMES else "0"
        start = str(row['start'])
        end = str(row['end'])
        geno = [re.split('\||/', x) for x in row['gts'].split(",")]
        geno = [["0", "0"] if any([y in self.NULL_GENOTYPES for y in x])
                else x for x in geno]
        genotypes = " ".join(list(flatten(geno)))
        alleles = "|".join(set(list(flatten(geno))).difference("0"))
        name = chrom + ":" +  start + "-" + end + ":" + alleles
        return " ".join([chrom, name, "0", start, genotypes])

    @classmethod
    def format_query(self, query):
        NEED_COLUMNS = ["chrom", "rs_ids", "start", "gts", "type", "variant_id"]
        return ensure_columns(query, NEED_COLUMNS)

    @classmethod
    def predicate(self, row):
        geno = [re.split("\||/", x) for x in row['gts']]
        geno = list(flatten(geno))
        num_alleles = len(set(geno).difference(self.NULL_GENOTYPES))
        return num_alleles > 0 and num_alleles <= 2 and row['type'] != "sv"


class JSONRowFormat(RowFormat):

    name = "json"

    @classmethod
    def format(self, row):
        """Emit a JSON representation of a given row
        """
        return json.dumps(row)

    @classmethod
    def format_query(self, query):
        return query

    @classmethod
    def predicate(self, row):
        return True


class GeminiRow(object):

    def __init__(self, row, gts=None, gt_types=None,
                 gt_phases=None, gt_depths=None,
                 gt_ref_depths=None, gt_alt_depths=None,
                 gt_quals=None, variant_samples=None,
                 HET_samples=None, HOM_ALT_samples=None,
                 formatter=DefaultRowFormat):
        self.row = row
        self.gts = gts
        self.gt_types = gt_types
        self.gt_phases = gt_phases
        self.gt_depths = gt_depths
        self.gt_ref_depths = gt_ref_depths
        self.gt_alt_depths = gt_alt_depths
        self.gt_quals = gt_quals
        self.gt_cols = ['gts', 'gt_types', 'gt_phases',
                        'gt_depths', 'gt_ref_depths', 'gt_alt_depths',
                        'gt_quals', "variant_samples", "HET_samples", "HOM_ALT_samples"]
        self.formatter = formatter
        self.variant_samples = variant_samples
        self.HET_samples = HET_samples
        self.HOM_ALT_samples = HOM_ALT_samples

    def __getitem__(self, val):
        if val not in self.gt_cols:
            return self.row[val]
        else:
            return getattr(self, val)

    def __iter__(self):
        return self

    def __repr__(self):
        return self.formatter.format(self.row)

    def next(self):
        try:
            return self.row.keys()
        except:
            raise StopIteration


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

    def __init__(self, db, include_gt_cols=False, out_format="default"):
        assert os.path.exists(db), "%s does not exist." % db

        self.db = db
        self.query_executed = False
        self.for_browser = False
        self.include_gt_cols = include_gt_cols

        self._connect_to_database()
        # map sample names to indices. e.g. self.sample_to_idx[NA20814] -> 323
        self.sample_to_idx = util.map_samples_to_indicies(self.c)
        # and vice versa. e.g., self.idx_to_sample[323] ->  NA20814
        self.idx_to_sample = util.map_indicies_to_samples(self.c)
        self.formatter = self._set_formatter(out_format.lower())
        self.predicates = [self.formatter.predicate]


    def _set_formatter(self, out_format):
        SUPPORTED_FORMATS = {x.name.lower(): x for x in itersubclasses(RowFormat)}

        if not out_format in SUPPORTED_FORMATS:
            raise NotImplementedError("Conversion to %s not supported. Valid "
                                      "formats are %s."
                                      % (out_format, SUPPORTED_FORMATS))
        else:
            return SUPPORTED_FORMATS[out_format]


    def _set_gemini_browser(self, for_browser):
        self.for_browser = for_browser

    def run(self, query, gt_filter=None, show_variant_samples=False,
            variant_samples_delim=',', predicates=None,
            needs_genotypes=False):
        """
        Execute a query against a Gemini database. The user may
        specify:

            1. (reqd.) an SQL `query`.
            2. (opt.) a genotype filter.
        """
        self.query = self.formatter.format_query(query)
        self.gt_filter = gt_filter
        self.show_variant_samples = show_variant_samples
        self.variant_samples_delim = variant_samples_delim
        if predicates:
            self.predicates += predicates

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
        if self.formatter.name in ["json"]:
            return None

        if self.query_type == "no-genotypes":
            h = [col for col in self.all_query_cols]
        else:
            h = [col for col in self.all_query_cols] + \
                [col for col in OrderedSet(self.all_columns_orig)
                 - OrderedSet(self.select_columns)]
        if self.show_variant_samples:
            h += ["variant_samples", "HET_samples", "HOM_ALT_samples"]
        return GeminiRow(OrderedDict(itertools.izip(h, h)),
                         formatter=self.formatter)

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
            except Exception as e:
                raise StopIteration
            gts = None
            gt_types = None
            gt_phases = None
            gt_depths = None
            gt_ref_depths = None
            gt_alt_depths = None
            gt_quals = None
            variant_names = []
            het_names = []
            hom_alt_names = []

            if self._query_needs_genotype_info():
                gts = compression.unpack_genotype_blob(row['gts'])
                gt_types = \
                    compression.unpack_genotype_blob(row['gt_types'])
                gt_phases = \
                    compression.unpack_genotype_blob(row['gt_phases'])
                gt_depths = \
                    compression.unpack_genotype_blob(row['gt_depths'])
                gt_ref_depths = \
                    compression.unpack_genotype_blob(row['gt_ref_depths'])
                gt_alt_depths = \
                    compression.unpack_genotype_blob(row['gt_alt_depths'])
                gt_quals = \
                    compression.unpack_genotype_blob(row['gt_quals'])
                variant_samples = [x for x, y in enumerate(gt_types) if y == HET or
                                   y == HOM_ALT]
                variant_names = [self.idx_to_sample[x] for x in variant_samples]
                het_samples = [x for x, y in enumerate(gt_types) if y == HET]
                het_names = [self.idx_to_sample[x] for x in het_samples]
                hom_alt_samples = [x for x, y in enumerate(gt_types) if y == HOM_ALT]
                hom_alt_names = [self.idx_to_sample[x] for x in hom_alt_samples]

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
                        elif col == "gt_quals":
                            fields[col] = \
                                ','.join(str(d) for d in gt_quals)
                        elif col == "gt_ref_depths":
                            fields[col] = \
                                ','.join(str(d) for d in gt_ref_depths)
                        elif col == "gt_alt_depths":
                            fields[col] = \
                                ','.join(str(d) for d in gt_alt_depths)

            if self.show_variant_samples:
                fields["variant_samples"] = \
                    self.variant_samples_delim.join(variant_names)
                fields["HET_samples"] = \
                    self.variant_samples_delim.join(het_names)
                fields["HOM_ALT_samples"] = \
                    self.variant_samples_delim.join(hom_alt_names)

            gemini_row = GeminiRow(fields, gts, gt_types, gt_phases,
                                   gt_depths, gt_ref_depths, gt_alt_depths,
                                   gt_quals, variant_names, het_names, hom_alt_names,
                                   formatter=self.formatter)

            if not all([predicate(gemini_row) for predicate in self.predicates]):
                continue

            if not self.for_browser:
                return gemini_row
            else:
                return fields

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

    def _execute_query(self):
        try:
            self.c.execute(self.query)
        except sqlite3.OperationalError as e:
            print "SQLite error: {0}".format(e)
            sys.exit("The query issued (%s) has a syntax error." % self.query)

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

            self._execute_query()

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
            self._execute_query()
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
            (column, sample) = raw_col.split('.', 1)
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

        if "from" not in self.query.lower():
            sys.exit("Malformed query: expected a FROM keyword.")

        (select_tokens, rest_of_query) = get_select_cols_and_rest(self.query)

        # remove any GT columns
        select_clause_list = []
        for token in select_tokens:
            if not token.startswith("gt") and not token.startswith("GT"):
                select_clause_list.append(token)

        # reconstruct the query with the GT* columns added
        if len(select_clause_list) > 0:
            select_clause = ",".join(select_clause_list) + \
                    ", gts, gt_types, gt_phases, gt_depths, \
                       gt_ref_depths, gt_alt_depths, gt_quals "

        else:
            select_clause = ",".join(select_clause_list) + \
                    " gts, gt_types, gt_phases, gt_depths, \
                      gt_ref_depths, gt_alt_depths, gt_quals "

        self.query = "select " + select_clause + rest_of_query

        # extract the original select columns
        return self.query


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
        if "from" not in self.query.lower():
            sys.exit("Malformed query: expected a FROM keyword.")

        (select_tokens, rest_of_query) = get_select_cols_and_rest(self.query)

        for token in select_tokens:
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
        return requested_genotype or \
               self.include_gt_cols or \
               self.show_variant_samples

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
