#!/usr/bin/env python

import os
import sys
import sqlite3
import re
import compiler
import collections
import json
import abc
import numpy as np
from itertools import chain
flatten = chain.from_iterable
import itertools as it

# gemini imports
import gemini_utils as util
from gemini_constants import *
from gemini_utils import (OrderedSet, OrderedDict, itersubclasses)
from .pdict import PDict
import compression
from sql_utils import ensure_columns, get_select_cols_and_rest
from gemini_subjects import get_subjects

class GeminiError(Exception):
    pass

def RowFactory(cursor, row):
    return dict(it.izip((c[0] for c in cursor.description), row))

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
        return str(row.print_fields)

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

    @abc.abstractmethod
    def header(self, fields):
        """ return a header for the row """
        return "\t".join(fields)


class DefaultRowFormat(RowFormat):

    name = "default"

    def __init__(self, args):
        pass

    def format(self, row):
        r = row.print_fields
        return '\t'.join(str(v) if not isinstance(v, np.ndarray) else ",".join(map(str, v)) for v in r._vals)

    def format_query(self, query):
        return query

    def predicate(self, row):
        return True

    def header(self, fields):
        """ return a header for the row """
        return "\t".join(fields)


class CarrierSummary(RowFormat):
    """
    Generates a count of the carrier/noncarrier status of each feature in a given
    column of the sample table

    Assumes None == unknown.
    """
    name = "carrier_summary"

    def __init__(self, args):
        subjects = get_subjects(args)
        self.carrier_summary = args.carrier_summary

        # get the list of all possible values in the column
        # but don't include None, since we are treating that as unknown.
        self.column_types = list(set([getattr(x, self.carrier_summary)
                                      for x in subjects.values()]))
        self.column_types = [i for i in self.column_types if i is not None]
        self.column_counters = {None: set()}
        for ct in self.column_types:
            self.column_counters[ct] = set([k for (k, v) in subjects.items() if
                                            getattr(v, self.carrier_summary) == ct])


    def format(self, row):
        have_variant = set(row['variant_samples'])
        have_reference = set(row['hom_ref_samples'])
        unknown = len(set(row['unknown_samples']).union(self.column_counters[None]))
        carrier_counts = []
        for ct in self.column_types:
            counts = len(self.column_counters[ct].intersection(have_variant))
            carrier_counts.append(counts)
        for ct in self.column_types:
            counts = len(self.column_counters[ct].intersection(have_reference))
            carrier_counts.append(counts)

        carrier_counts.append(unknown)
        carrier_counts = map(str, carrier_counts)
        r = row.print_fields
        return '\t'.join(str(r[c]) if not isinstance(r[c], np.ndarray) else
                         ",".join(str(s) for s in r[c]) for c in r) \
               + "\t" + "\t".join(carrier_counts)

    def format_query(self, query):
        return query

    def predicate(self, row):
        return True

    def header(self, fields):
        """ return a header for the row """
        header_columns = self.column_types
        if self.carrier_summary == "affected":
            header_columns = self._rename_affected()
        carriers = [x + "_carrier" for x in map(str, header_columns)]
        noncarriers = [x + "_noncarrier" for x in map(str, header_columns)]
        fields += carriers
        fields += noncarriers
        fields += ["unknown"]
        return "\t".join(fields)

    def _rename_affected(self):
        header_columns = []
        for ct in self.column_types:
            if ct is True:
                header_columns.append("affected")
            elif ct is False:
                header_columns.append("unaffected")
        return header_columns



class TPEDRowFormat(RowFormat):

    X_PAR_REGIONS = [(60001, 2699520), (154931044, 155260560)]
    Y_PAR_REGIONS = [(10001, 2649520), (59034050, 59363566)]

    name = "tped"
    NULL_GENOTYPES = ["."]
    PED_MISSING = ["0", "0"]
    VALID_CHROMOSOMES = map(str, range(1, 23)) + ["X", "Y", "XY", "MT"]
    POSSIBLE_HAPLOID = ["X", "Y"]

    def __init__(self, args):
        gq = GeminiQuery(args.db)
        subjects = get_subjects(args, skip_filter=True)
        # get samples in order of genotypes
        self.samples = [gq.idx_to_sample_object[x] for x in range(len(subjects))]

    def format(self, row):
        VALID_CHROMOSOMES = map(str, range(1, 23)) + ["X", "Y", "XY", "MT"]
        chrom = row['chrom'].split("chr")[1]
        chrom = chrom if chrom in VALID_CHROMOSOMES else "0"
        start = str(row.row['start'])
        end = str(row.row['end'])
        ref = row['ref']
        alt = row['alt']
        geno = [re.split('\||/', x) for x in row['gts']]
        geno = [self._fix_genotype(chrom, start, genotype, self.samples[i].sex)
                for i, genotype in enumerate(geno)]
        genotypes = " ".join(list(flatten(geno)))
        name = str(row['variant_id'])
        return " ".join([chrom, name, "0", start, genotypes])

    def format_query(self, query):
        NEED_COLUMNS = ["chrom", "rs_ids", "start", "ref", "alt", "gts", "type", "variant_id"]
        return ensure_columns(query, NEED_COLUMNS)


    def predicate(self, row, _splitter=re.compile("\||/")):
        geno = [_splitter.split(x) for x in row['gts']]
        geno = list(flatten(geno))
        num_alleles = len(set(geno).difference(self.NULL_GENOTYPES))
        return num_alleles > 0 and num_alleles <= 2 and row['type'] != "sv"

    def _is_haploid(self, genotype):
        return len(genotype) < 2

    def _has_missing(self, genotype):
        return any([allele in self.NULL_GENOTYPES for allele in genotype])

    def _is_heterozygote(self, genotype):
        return len(genotype) == 2 and (genotype[0] != genotype[1])

    def _in_PAR(self, chrom, start):
        if chrom == "X":
            for region in self.X_PAR_REGIONS:
                if start > region[0] and start < region[1]:
                    return True
        elif chrom == "Y":
            for region in self.Y_PAR_REGIONS:
                if start > region[0] and start < region[1]:
                    return True
        return False

    def _fix_genotype(self, chrom, start, genotype, sex):
        """
        the TPED format has to have both alleles set, even if it is haploid.
        this fixes that setting Y calls on the female to missing,
        heterozygotic calls on the male non PAR regions to missing and haploid
        calls on non-PAR regions to be the haploid call for both alleles
        """
        if sex == "2":
            # set female Y calls and haploid calls to missing
            if self._is_haploid(genotype) or chrom == "Y" or self._has_missing(genotype):
                return self.PED_MISSING
            return genotype
        if chrom in self.POSSIBLE_HAPLOID and sex == "1":
            # remove the missing genotype calls
            genotype = [x for x in genotype if x not in self.NULL_GENOTYPES]
            # if all genotypes are missing skip
            if not genotype:
                return self.PED_MISSING
            # heterozygote males in non PAR regions are a mistake
            if self._is_heterozygote(genotype) and not self._in_PAR(chrom, start):
                return self.PED_MISSING
            # set haploid males to be homozygous for the allele
            if self._is_haploid(genotype):
                return [genotype[0], genotype[0]]

        # if a genotype is missing or is haploid set it to missing
        if self._has_missing(genotype) or self._is_haploid(genotype):
            return self.PED_MISSING
        else:
            return genotype

    def header(self, fields):
        return None

class JSONRowFormat(RowFormat):

    name = "json"

    def __init__(self, args):
        pass

    def format(self, row):
        """Emit a JSON representation of a given row
        """
        from .pdict import to_json
        return json.dumps(row.print_fields, default=to_json)

    def format_query(self, query):
        return query

    def predicate(self, row):
        return True

    def header(self, fields):
        return None

class VCFRowFormat(RowFormat):

    name = "vcf"

    def __init__(self, args):
        self.gq = GeminiQuery(args.db)

    def format(self, row):
        """Emit a VCF representation of a given row

           TODO: handle multiple alleles
        """
        # core VCF fields
        vcf_rec = [row.row['chrom'], row.row['start'] + 1]
        if row.row['vcf_id'] is None:
            vcf_rec.append('.')
        else:
            vcf_rec.append(row.row['vcf_id'])
        vcf_rec += [row.row['ref'], row.row['alt'], row.row['qual']]
        if row.row['filter'] is None:
            vcf_rec.append('PASS')
        else:
            vcf_rec.append(row.row['filter'])
        vcf_rec += [row.row['info'], 'GT']

        # construct genotypes
        gts = list(row['gts'])
        gt_types = list(row['gt_types'])
        gt_phases = list(row['gt_phases'])
        for idx, gt_type in enumerate(gt_types):
            phase_char = '/' if not gt_phases[idx] else '|'
            gt = gts[idx]
            alleles = gt.split(phase_char)
            if gt_type == HOM_REF:
                vcf_rec.append('0' + phase_char + '0')
            elif gt_type == HET:
                # if the genotype is phased, need to check for 1|0 vs. 0|1
                if gt_phases[idx] and alleles[0] != row.row['ref']:
                    vcf_rec.append('1' + phase_char + '0')
                else:
                    vcf_rec.append('0' + phase_char + '1')
            elif gt_type == HOM_ALT:
                vcf_rec.append('1' + phase_char + '1')
            elif gt_type == UNKNOWN:
                vcf_rec.append('.' + phase_char + '.')

        return '\t'.join([str(c) if c is not None else "." for c in vcf_rec])

    def format_query(self, query):
        return query

    def predicate(self, row):
        return True

    def header(self, fields):
        """Return the original VCF's header
        """
        try:
            self.gq.run('select vcf_header from vcf_header')
            return str(self.gq.next()).strip()
        except:
            sys.exit("Your database does not contain the vcf_header table. Therefore, you cannot use --header.\n")

class SampleDetailRowFormat(RowFormat):
    """Retrieve queries with flattened sample information for samples present.

    This melts/tidys a single line result to have a separate line for every sample
    with that call and adds in sample metadata from the samples table.
    """
    name = "sampledetail"

    def __init__(self, args):
        self.gq = GeminiQuery(args.db)
        self.gq.run("SELECT * from samples")
        self.cols = self.gq.header.split()[1:]
        self.args = args

        self.samples = {}
        for row in self.gq:
            vals = [row[x] for x in self.cols]
            self.samples[row["name"]] = vals

    def format(self, row):
        r = row.print_fields
        samples = [s for s in r["variant_samples"].split(self.args.sample_delim) if s]
        for x in self.gq.sample_show_fields:
            del r[x]
        out = []
        for sample in samples:
            out.append('\t'.join([str(r[c]) for c in r] + self.samples[sample]))
        return "\n".join(out)

    def format_query(self, query):
        return query

    def predicate(self, row):
        return True

    def header(self, fields):
        for x in self.gq.sample_show_fields:
            if x in fields:
                fields.remove(x)
        return "\t".join(fields + self.cols)

class GeminiRow(object):
    __slots__ = ('cache', 'genotype_dict', 'row', 'formatter', 'query',
                 'print_fields')

    def __init__(self, row, query, formatter=DefaultRowFormat(None), print_fields=None):
        # row can be a dict() from the database or another GeminiRow (from the
        # same db entry). we try to re-use the cached stuff if possible.
        self.cache = {}
        self.genotype_dict = {}
        self.row = getattr(row, "row", row)
        self.cache = getattr(row, "cache", {})
        self.genotype_dict = getattr(row, "genotype_dict", {})

        # for the eval.
        #self.cache['sample_info'] = dict(query.sample_info)

        self.formatter = formatter
        self.query = query
        self.print_fields = print_fields or {}

    def __getitem__(self, key):
        # we cache what we can.
        key = str(key)
        if key in ('het_samples', 'hom_alt_samples', 'unknown_samples',
                'variant_samples', 'hom_ref_samples'):
            if self.genotype_dict == {}:
                self.genotype_dict = self.query._group_samples_by_genotype(self['gt_types'])
            if key == 'het_samples':
                return self.genotype_dict[HET]
            if key == 'hom_alt_samples':
                return self.genotype_dict[HOM_ALT]
            if key == 'hom_ref_samples':
                return self.genotype_dict[HOM_REF]
            if key == 'unknown_samples':
                return self.genotype_dict[UNKNOWN]
            if key == 'variant_samples':
                return self.genotype_dict[HET] + self.genotype_dict[HOM_ALT]

        if key in self.cache:
            return self.cache[key]

        if key == 'info':
            if 'info' not in self.cache:
                self.cache['info'] = compression.unpack_ordereddict_blob(self.row['info'])
            return self.cache['info']
        if key not in self.query.gt_cols:
            return self.row[key]
        elif key in self.query.gt_cols:
            if key not in self.cache:
                self.cache[key] = compression.unpack_genotype_blob(self.row[key])
            return self.cache[key]
        raise KeyError(key)

    def keys(self):
        return self.row.keys() + self.cache.keys()

    def __iter__(self):
        return self

    def __repr__(self):
        return self.formatter.format(self)

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

    def __init__(self, db, include_gt_cols=False,
                 out_format=DefaultRowFormat(None),
                 variant_id_getter=None):
        assert os.path.exists(db), "%s does not exist." % db

        self.db = db
        self.query_executed = False
        self.for_browser = False
        self.include_gt_cols = include_gt_cols
        self.variant_id_getter = variant_id_getter

        # try to connect to the provided database
        self._connect_to_database()

        # save the gt_cols in the database and don't hard-code them anywhere.
        self.gt_cols = util.get_gt_cols(self.conn)

        # extract the column names from the sample table.
        # needed for gt-filter wildcard support.
        self._collect_sample_table_columns()

        # list of samples ids for each clause in the --gt-filter
        self.sample_info = collections.defaultdict(list)

        # map sample names to indices. e.g. self.sample_to_idx[NA20814] -> 323
        self.sample_to_idx = util.map_samples_to_indices(self.c)
        # and vice versa. e.g., self.idx_to_sample[323] ->  NA20814
        self.idx_to_sample = util.map_indices_to_samples(self.c)
        self.idx_to_sample_object = util.map_indices_to_sample_objects(self.c)
        self.sample_to_sample_object = util.map_samples_to_sample_objects(self.c)
        self.formatter = out_format
        self.predicates = [self.formatter.predicate]
        self.sample_show_fields = ["variant_samples", "het_samples", "hom_alt_samples"]

    def _set_gemini_browser(self, for_browser):
        self.for_browser = for_browser

    def run(self, query, gt_filter=None, show_variant_samples=False,
            variant_samples_delim=',', predicates=None,
            needs_genotypes=False, needs_genes=False,
            show_families=False, subjects=None):
        """
        Execute a query against a Gemini database. The user may
        specify:

            1. (reqd.) an SQL `query`.
            2. (opt.) a genotype filter.
        """
        self.query = self.formatter.format_query(query)
        self.gt_filter = gt_filter
        if self._is_gt_filter_safe() is False:
            sys.exit("ERROR: invalid --gt-filter command.")

        self.show_variant_samples = show_variant_samples
        self.variant_samples_delim = variant_samples_delim

        self.needs_genotypes = needs_genotypes
        self.needs_vcf_columns = False
        if self.formatter.name == 'vcf':
            self.needs_vcf_columns = True

        self.needs_genes = needs_genes
        self.show_families = show_families
        self.subjects = subjects
        if predicates:
            self.predicates += predicates

        # make sure the SELECT columns are separated by a
        # comma and a space. then tokenize by spaces.
        self.query = self.query.replace(',', ', ')
        self.query_pieces = self.query.split()
        if not any(s.startswith(("gt", "(gt")) for s in self.query_pieces) and \
           not any(".gt" in s for s in self.query_pieces):
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

        if self.gt_filter:
            # here's how we use the fast
            if self.variant_id_getter:
                if os.environ.get('GEMINI_DEBUG') == 'TRUE':
                    print >>sys.stderr, "bcolz: using index"

                user_dict = dict(HOM_REF=0, HET=1, UNKNOWN=2, HOM_ALT=3,
                                 sample_info=self.sample_info,
                                 MISSING=None, UNAFFECTED=1, AFFECTED=2)
                import time
                t0 = time.time()
                vids = self.variant_id_getter(self.db, self.gt_filter, user_dict)
                if vids is None:
                    print >>sys.stderr, "bcolz: can't parse this filter (falling back to gemini): %s" % self.gt_filter
                else:
                    if os.environ.get('GEMINI_DEBUG') == 'TRUE':
                        print >>sys.stderr, "bcolz: %.2f seconds to get %d rows." % (time.time() - t0, len(vids))
                    self.add_vids_to_query(vids)

        if self.gt_filter:
            self.gt_filter_compiled = compiler.compile(self.gt_filter, self.gt_filter, 'eval')

        self._apply_query()
        self.query_executed = True

    def __iter__(self):
        return self

    def add_vids_to_query(self, vids):
        #extra = " variant_id IN (%s)" % ",".join(map(str, vids))
        # faster way to convert to string.
        # http://stackoverflow.com/a/13861407
        q = add_variant_ids_to_query(self.query, vids)
        if q:
            self.query = q
            self.gt_filter = None

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
            h += self.sample_show_fields
        if self.show_families:
            h += ["families"]
        return self.formatter.header(h)

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
                row = GeminiRow(self.c.next(), self)
            except Exception:
                self.conn.close()
                raise StopIteration

            # skip the record if it does not meet the user's genotype filter
            # short circuit some expensive ops
            if self.gt_filter:
                try:
                    if 'False' == self.gt_filter: continue
                    unpacked = {'sample_info': self.sample_info}
                    for col in self.gt_cols:
                        if col in self.gt_filter:
                            unpacked[col] = row[col]

                    if not eval(self.gt_filter_compiled, unpacked):
                        continue
                # eval on a phred_ll column that was None
                except TypeError:
                    continue

            fields = PDict()

            for idx, col in enumerate(self.report_cols):
                if col == "*":
                    continue
                if not col[:2] in ("gt", "GT"):
                    # need to use add in case of duplicated fields (from
                    # variants and variant_impacts)
                    fields.add(col, row[col])
                else:
                    # reuse the original column name user requested
                    # e.g. replace gts[1085] with gts.NA20814
                    if '[' in col:
                        orig_col = self.gt_idx_to_name_map[col]

                        source, extra = col.split('[', 1)
                        assert extra[-1] == ']'
                        if source.startswith('gt_phred_ll') and row[source] is None:
                            fields[orig_col] = None
                            continue

                        idx = int(extra[:-1])
                        val = row[source][idx]

                        if type(val) in (np.int8, np.int32, np.bool_):
                            fields.add(orig_col, int(val))

                        elif type(val) in (np.float32,):
                            fields.add(orig_col, float(val))
                        else:
                            fields.add(orig_col, val)
                    else:
                        # asked for "gts" or "gt_types", e.g.
                        if row[col] is not None:
                            fields[col] = row[col]
                        else:
                            fields[col] = str(None)

            if self.show_variant_samples:
                fields["variant_samples"] = \
                    self.variant_samples_delim.join(self._filter_samples(row['variant_samples']))
                fields["het_samples"] = \
                    self.variant_samples_delim.join(self._filter_samples(row['het_samples']))
                fields["hom_alt_samples"] = \
                    self.variant_samples_delim.join(self._filter_samples(row['hom_alt_samples']))
            if self.show_families:
                families = map(str, list(set([self.sample_to_sample_object[x].family_id
                                              for x in row['variant_samples']])))
                fields["families"] = self.variant_samples_delim.join(families)

            if not all(predicate(row) for predicate in self.predicates):
                continue

            if not self.for_browser:
                # need to use new row for formatter.
                return GeminiRow(row, self, formatter=self.formatter, print_fields=fields)
            else:
                return fields

    def _filter_samples(self, samples):
        """Respect --sample-filter when outputting lists of sample information.
        """
        if self.subjects is not None:
            return [x for x in samples if x in self.subjects]
        else:
            return samples

    def _group_samples_by_genotype(self, gt_types):
        """
        make list keyed by genotype of list of samples with that genotype
        so index 0 is HOM, 1 is HET, 2 is UKNOWN, 3 is HOM_ALT.
        """
        d = [[], [], [], []]
        lookup = self.idx_to_sample
        for i, x in enumerate(gt_types):
            d[x].append(lookup[i])
        return d

    def _connect_to_database(self):
        """
        Establish a connection to the requested Gemini database.
        """
        # open up a new database
        if os.path.exists(self.db):
            self.conn = sqlite3.connect(self.db)
            self.conn.isolation_level = None
            # allow us to refer to columns by name
            #self.conn.row_factory = RowFactory
            self.conn.row_factory = sqlite3.Row
            self.c = self.conn.cursor()


    def _collect_sample_table_columns(self):
        """
        extract the column names in the samples table into a list
        """
        self.c.execute('select * from samples limit 1')
        self.sample_column_names = [tup[0] for tup in self.c.description]

    def _is_gt_filter_safe(self):
        """
        Test to see if the gt_filter string is potentially malicious.

        A future improvement would be to use pyparsing to
        traverse and directly validate the string.
        """
        if self.gt_filter is None or len(self.gt_filter.strip()) == 0:
            return True

        # avoid builtins
        # http://nedbatchelder.com/blog/201206/eval_really_is_dangerous.html
        if "__" in self.gt_filter:
            return False

        # avoid malicious commands
        evil = [" rm ", "os.system"]
        if any(s in self.gt_filter for s in evil):
            return False

        # make sure a "gt" col is in the string
        valid_cols = list(flatten(("%s." % gtc, "(%s)." % gtc) for gtc in self.gt_cols))
        if any(s in self.gt_filter for s in valid_cols):
            return True

        # assume the worst
        return False

    def _execute_query(self):
        try:
            self.c.execute(self.query)
        except sqlite3.OperationalError as e:
            msg = "SQLite error: {0}\n".format(e)
            print msg
            sys.stderr.write(msg)
            sys.exit("The query issued (%s) has a syntax error." % self.query)

    def _apply_query(self):
        """
        Execute a query. Intercept gt* columns and
        replace sample names with indices where necessary.
        """
        if self.needs_genes:
            self.query = self._add_gene_col_to_query()

        if self.needs_vcf_columns:
            self.query = self._add_vcf_cols_to_query()

        if self._query_needs_genotype_info():
            # break up the select statement into individual
            # pieces and replace genotype columns using sample
            # names with sample indices
            self._split_select()

            # we only need genotype information if the user is
            # querying the variants table
            self.query = self._add_gt_cols_to_query()
            self._execute_query()

            self.all_query_cols = [
                str(tuple[0]) for tuple in self.c.description
                if not tuple[0][:2] == "gt" and ".gt" not in tuple[0]
                ]

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
                    if not tuple[0][:2] == "gt"]
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
            return raw_col
        # e.g., "gts.NA12878"
        elif '.' in raw_col:
            (column, sample) = raw_col.split('.', 1)
            corrected = "%s[%d]" % (column.lower(), self.sample_to_idx[sample])
        else:
            # e.g. "gts" - do nothing
            corrected = raw_col
        return corrected

    def _get_matching_sample_ids(self, wildcard):
        """
        Helper function to convert a sample wildcard
        to a list of tuples reflecting the sample indices
        and sample names so that the wildcard
        query can be applied to the gt_* columns.
        """
        query = 'SELECT sample_id, name FROM samples '
        if wildcard.strip() != "*":
           query += ' WHERE ' + wildcard

        sample_info = [] # list of sample_id/name tuples
        self.c.execute(query)
        for row in self.c:
            # sample_ids are 1-based but gt_* indices are 0-based
            sample_info.append((int(row['sample_id']) - 1, str(row['name'])))
        return sample_info

    def _correct_genotype_filter(self):
        """
        This converts a raw genotype filter that contains
        'wildcard' statements into a filter that can be eval()'ed.
        Specifically, we must convert a _named_ genotype index
        to a _numerical_ genotype index so that the appropriate
        value can be extracted for the sample from the genotype
        numpy arrays.

        For example, without WILDCARDS, this converts:
        --gt-filter "(gt_types.1478PC0011 == 1)"

        to:
        (gt_types[11] == 1)

        With WILDCARDS, this converts things like:
            "(gt_types).(phenotype==1).(==HET)"

        to:
            "gt_types[2] == HET and gt_types[5] == HET"
        """

        def _swap_genotype_for_number(token):
            """
            This is a bit of a hack to get around the fact that eval()
            doesn't handle the imported constants well when also having to
            find local variables.  This requires some eval/globals()/locals() fu
            that has evaded me thus far. Just replacing HET, etc. with 1, etc. works.
            """
            if any(g in token for g in ['HET', 'HOM_ALT', 'HOM_REF', 'UNKNOWN']):
                token = token.replace('HET', str(HET))
                token = token.replace('HOM_ALT', str(HOM_ALT))
                token = token.replace('HOM_REF', str(HOM_REF))
                token = token.replace('UNKNOWN', str(UNKNOWN))
            return token

        corrected_gt_filter = []

        # first try to identify wildcard rules.
        # (\s*gt\w+\) handles both
        #    (gt_types).(*).(!=HOM_REF).(all)
        # and
        #    (   gt_types).(*).(!=HOM_REF).(all)
        seen_count = False
        wildcard_tokens = re.split(r'(\(\s*gt\w+\s*\)\.\(.+?\)\.\(.+?\)\.\(.+?\))', str(self.gt_filter))
        for token_idx, token in enumerate(wildcard_tokens):
            # NOT a WILDCARD
            # We must then split on whitespace and
            # correct the gt_* columns:
            # e.g., "gts.NA12878" or "and gt_types.M10500 == HET"
            if (token.find("gt") >= 0 or token.find("GT") >= 0) \
                and not '.(' in token and not ')self.' in token:
                tokens = re.split(r'[\s+]+', str(token))
                for t in tokens:
                    if len(t) == 0:
                        continue
                    if (t.find("gt") >= 0 or t.find("GT") >= 0):
                        corrected = self._correct_genotype_col(t)
                        corrected_gt_filter.append(corrected)
                    else:
                        t = _swap_genotype_for_number(t)
                        if t.strip() in ("AND", "OR"):
                            t = t.lower()
                        corrected_gt_filter.append(t)
            # IS a WILDCARD
            # e.g., "gt_types.(affected==1).(==HET)"
            elif (token.find("gt") >= 0 or token.find("GT") >= 0) \
                and '.(' in token and ').' in token:
                # break the wildcard into its pieces. That is:
                # (COLUMN).(WILDCARD).(WILDCARD_RULE).(WILDCARD_OP)
                # e.g, (gts).(phenotype==2).(==HET).(any)
                if token.count('.') != 3 or \
                   token.count('(') != 4 or \
                   token.count(')') != 4:
                    sys.exit("Wildcard filter should consist of 4 elements. Exiting.")

                (column, wildcard, wildcard_rule, wildcard_op) = token.split('.')

                # remove the syntactic parentheses
                column = column.strip('(').strip(')').strip()
                wildcard = wildcard.strip('(').strip(')').strip()
                wildcard_rule = wildcard_rule.strip('(').strip(')').strip()
                wildcard_op = wildcard_op.strip('(').strip(')').strip()

                # collect and save all of the samples that meet the wildcard criteria
                # for each clause.
                # these will be used in the list comprehension for the eval expression
                # constructed below.
                self.sample_info[token_idx] = self._get_matching_sample_ids(wildcard)

                # Replace HET, etc. with 1, et.c to avoid eval() issues.
                wildcard_rule = _swap_genotype_for_number(wildcard_rule)

                # build the rule based on the wildcard the user has supplied.
                if wildcard_op in ["all", "any"]:
                    if self.variant_id_getter:
                        joiner = " and " if wildcard_op == "all" else " or "
                        rule = joiner.join("%s[%s]%s" % (column, s[0], wildcard_rule) for s in self.sample_info[token_idx])
                        rule = "(" + rule + ")"
                    else:
                        rule = wildcard_op + "(" + column + '[sample[0]]' + wildcard_rule + " for sample in sample_info[" + str(token_idx) + "])"
                elif wildcard_op == "none":
                    if self.variant_id_getter:
                        rule = " or ".join("%s[%s]%s" % (column, s[0], wildcard_rule) for s in self.sample_info[token_idx])
                        rule = "~ ((" + rule + "))"
                    else:
                        rule = "not any(" + column + '[sample[0]]' + wildcard_rule + " for sample in sample_info[" + str(token_idx) + "])"
                elif "count" in wildcard_op:
                    # break "count>=2" into ['', '>=2']
                    tokens = wildcard_op.split('count')
                    count_comp = tokens[len(tokens) - 1]
                    if self.variant_id_getter:
                        rule = "|count|".join("((%s[%s]%s))" % (column, s[0], wildcard_rule) for s in self.sample_info[token_idx])
                        rule = "%s|count|%s" % (rule, count_comp.strip())
                        seen_count = True
                    else:
                        rule = "sum(" + column + '[sample[0]]' + wildcard_rule + " for sample in sample_info[" + str(token_idx) + "])" + count_comp
                else:
                    sys.exit("Unsupported wildcard operation: (%s). Exiting." % wildcard_op)

                corrected_gt_filter.append(rule)
            else:
                if len(token) > 0:
                    if token.strip() in ("AND", "OR"):
                        token = token.lower()
                    corrected_gt_filter.append(token)
        if seen_count and len(corrected_gt_filter) > 1 and self.variant_id_getter:
            raise GeminiError("count operations can not be combined with other operations")

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
            if not token[:2] in ("GT", "gt") and \
               not token[:3] in ("(gt", "(GT") and \
               not ".gt" in token and \
               not ".GT" in token:
                select_clause_list.append(token)

        # reconstruct the query with the GT* columns added
        select_clause = (", ".join(select_clause_list)).strip()
        if select_clause: select_clause += ","
        select_clause += ",".join(self.gt_cols)

        self.query = "select %s %s" % (select_clause, rest_of_query)

        # extract the original select columns
        return self.query

    def _add_gene_col_to_query(self):
        """
        Add the gene column to the list of SELECT'ed columns
        in a query.
        """
        if "from" not in self.query.lower():
            sys.exit("Malformed query: expected a FROM keyword.")

        (select_tokens, rest_of_query) = get_select_cols_and_rest(self.query)

        if not any("gene" in s for s in select_tokens):

            select_clause = ",".join(select_tokens) + \
                        ", gene "

            self.query = "select " + select_clause + rest_of_query

        return self.query

    def _add_vcf_cols_to_query(self):
        """
        Add the VCF columns to the list of SELECT'ed columns
        in a query.

        NOTE: Should only be called if using VCFRowFormat()
        """
        if "from" not in self.query.lower():
            sys.exit("Malformed query: expected a FROM keyword.")

        (select_tokens, rest_of_query) = get_select_cols_and_rest(self.query)

        cols_to_add = []
        for col in ['chrom', 'start', 'vcf_id', 'ref', 'alt', 'qual', 'filter', 'info', \
            'gts', 'gt_types', 'gt_phases']:
            if not any(col in s for s in select_tokens):
                cols_to_add.append(col)

        select_clause = ",".join(select_tokens + cols_to_add)
        self.query = "select " + select_clause + rest_of_query
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

            # it is a WILDCARD
            if (token.find("gt") >= 0 or token.find("GT") >= 0) \
                and '.(' in token and ').' in token:
                # break the wildcard into its pieces. That is:
                # (COLUMN).(WILDCARD)
                (column, wildcard) = token.split('.')

                # remove the syntactic parentheses
                wildcard = wildcard.strip('(').strip(')')
                column = column.strip('(').strip(')')

                # convert "gt_types.(affected==1)"
                # to: gt_types[3] == HET and gt_types[9] == HET
                sample_info = self._get_matching_sample_ids(wildcard)

                # maintain a list of the sample indices that should
                # be displayed as a result of the SELECT'ed wildcard
                wildcard_indices = []
                for (idx, sample) in enumerate(sample_info):
                    wildcard_display_col = column + '.' + str(sample[1])
                    wildcard_mask_col = column + '[' + str(sample[0]) + ']'
                    wildcard_indices.append(sample[0])

                    new_col = wildcard_mask_col
                    self.all_columns_new.append(new_col)
                    self.all_columns_orig.append(wildcard_display_col)
                    self.gt_name_to_idx_map[wildcard_display_col] = wildcard_mask_col
                    self.gt_idx_to_name_map[wildcard_mask_col] = wildcard_display_col

            # it is a basic genotype column
            elif (token.find("gt") >= 0 or token.find("GT") >= 0) \
                and '.(' not in token and not ').' in token \
                and "length" not in token:     # e.g., no aa_lenGTh, etc. false positives
                new_col = self._correct_genotype_col(token)

                self.all_columns_new.append(new_col)
                self.all_columns_orig.append(token)
                self.gt_name_to_idx_map[token] = new_col
                self.gt_idx_to_name_map[new_col] = token

            # it is neither
            else:
                self.select_columns.append(token)
                self.all_columns_new.append(token)
                self.all_columns_orig.append(token)

    def _info_dict_to_string(self, info):
        """
        Flatten the VCF info-field OrderedDict into a string,
        including all arrays for allelic-level info.
        """
        if info is not None:
            return ';'.join('%s=%s' % (key, value) if not isinstance(value, list)
                             else '%s=%s' % (key, ','.join(str(v) for v in value))
                             for (key, value) in info.items())
        else:
            return None

    def _tokenize_query(self):
        return list(flatten(x.split(",") for x in self.query.split(" ")))

    def _query_needs_genotype_info(self):
        if self.include_gt_cols or self.show_variant_samples or self.needs_genotypes:
            return True

        tokens = self._tokenize_query()
        requested_genotype = "variants" in tokens and \
                            (any(x.startswith(("gt", "(gt")) for x in tokens) or \
                             any(".gt" in x for x in tokens))
        return requested_genotype

def select_formatter(args):
    SUPPORTED_FORMATS = {x.name.lower(): x for x in
                         itersubclasses(RowFormat)}

    if hasattr(args, 'carrier_summary') and args.carrier_summary:
        return SUPPORTED_FORMATS["carrier_summary"](args)

    if args.format not in SUPPORTED_FORMATS:
        raise NotImplementedError("Conversion to %s not supported. Valid "
                                  "formats are %s."
                                  % (args.format, SUPPORTED_FORMATS))
    else:
        return SUPPORTED_FORMATS[args.format](args)

def add_variant_ids_to_query(query, vids):
    """
    >>> vids = range(1, 4)
    >>> add_variant_ids_to_query("select gene, chrom, start, end from variants limit 10", vids)
    'select gene, chrom, start, end from variants where  variant_id IN (1,2,3)  limit 10'
    >>> add_variant_ids_to_query("select gene, chrom, start, end from variants where gene = 'asdf' limit 10", vids)
    "select gene, chrom, start, end from variants where gene = 'asdf' and  variant_id IN (1,2,3)  limit 10"
    >>> add_variant_ids_to_query("select gene, chrom, start, end from variants where gene = 'asdf' order by gene limit 10", vids)
    "select gene, chrom, start, end from variants where gene = 'asdf' and  variant_id IN (1,2,3) order by gene limit 10"
    >>> add_variant_ids_to_query("select gene, chrom, start, end from variants", vids)
    'select gene, chrom, start, end from variants where  variant_id IN (1,2,3)'
    """
    if len(vids) == 0:
        return vids
    extra = " variant_id IN (%s)" % ",".join(np.char.mod("%i", vids))

    # order by, then limit.
    limit_idx = query.lower().index(" limit ") if " limit " in query.lower() else None
    if limit_idx:
        query, qlimit = query[:limit_idx].strip(), query[limit_idx:].strip()
        assert qlimit.lower().startswith("limit ")
    else:
        qlimit = ""

    order_idx = query.lower().index(" order by ") if " order by " in query.lower() else None
    if order_idx:
        query, qorder = query[:order_idx].strip(), query[order_idx:].strip()
        assert qorder.lower().startswith("order by ")
    else:
        qorder = ""

    if " where " in query.lower():
        extra = " and " + extra
    else:
        extra = " where " + extra
    return " ".join([x.strip() for x in (query, extra, qorder, qlimit)]).strip()


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
