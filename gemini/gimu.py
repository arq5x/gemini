#!/usr/bin/env python
from collections import Counter
import GeminiQuery
import sql_utils
import compiler
from gemini_constants import *
from .gemini_bcolz import filter
import itertools as it
import operator as op


class GeminiInheritanceModel(object):

    required_columns = ("family_id", "family_members",
                        "family_genotypes", "samples", "family_count")

    def __init__(self, args):

        self.args = args
        self.gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)

        self.gt_cols = self.gq.gt_cols

        if not args.columns:
            args.columns = "*," + ", ".join(self.gt_cols)
        self.set_family_info()

    @property
    def query(self):
        # TODO: need
        if self.args.columns is not None:
            # the user only wants to report a subset of the columns
            query = "SELECT " + str(self.args.columns) + " FROM variants "
        else:
            # report the kitchen sink
            query = "SELECT chrom, start, end, * %s " \
                + "FROM variants" % ", ".join(self.gt_cols)

        query = sql_utils.ensure_columns(query, ['variant_id'])
        # add any non-genotype column limits to the where clause
        if self.args.filter:
            query += " WHERE " + self.args.filter

        # auto_rec and auto_dom candidates should be limited to
        # variants affecting genes.
        if self.model in ("auto_rec", "auto_dom") or \
           (self.model == "de_novo" and self.args.min_kindreds is not None):

            # we require the "gene" column for the auto_* tools
            query = sql_utils.ensure_columns(query, ['gene'])
            if self.args.filter:
                query += " AND gene is not NULL ORDER BY chrom, gene"
            else:
                query += " WHERE gene is not NULL ORDER BY chrom, gene"
        return query

    def bcolz_candidates(self):
        """
        Get all the variant ids that meet the genotype filter for any fam.
        """
        variant_ids = set()
        for i, family_id in enumerate(self.family_ids):
            gt_filter = self.family_masks[i]
            variant_ids.update(filter(self.args.db, gt_filter, {}))

        return sorted(set(variant_ids))

    def gen_candidates(self, group_key):
        if isinstance(group_key, basestring):
            group_key = op.itemgetter(group_key)

        q = self.query
        vids = self.bcolz_candidates()
        q = GeminiQuery.add_variant_ids_to_query(q, vids)
        self.gq.run(q, needs_genotypes=True)

        def update(gr):
            # gr is a gemini row
            return gr

        for grp_key, grp in it.groupby(self.gq, group_key):
            ogrp = (update(gr) for gr in grp)
            yield grp_key, ogrp

    def all_candidates(self):

        _, candidates = self.gen_candidates(group_key=None)
        for candidate in candidates:
            yield candidate

    def gene_candidates(self):

        for gene, candidates in self.gen_candidates(group_key="gene"):
            yield gene, candidates

    def set_family_info(self):
        """
        Extract the relevant genotype filters, as well all labels
        for each family in the database.
        """
        from .family import Family
        self.families = families = Family.from_cursor(self.gq.c).values()

        self.family_ids = []
        self.family_masks = []
        for family in families:
            # e.g. family.auto_rec(gt_ll, min_depth)
            family_filter = getattr(family,
                    self.model)(gt_ll=self.args.gt_phred_ll,
                                min_depth=self.args.min_sample_depth)

            self.family_masks.append(family_filter)
            self.family_ids.append(family.family_id)

    def report_candidates(self):
        req_cols = ['gt_types', 'gts']
        if self.args.min_sample_depth and self.args.min_sample_depth > 0:
            req_cols.append('gt_depths')
        if self.args.gt_phred_ll:
            req_cols.extend(['gt_phred_ll_homref', 'gt_phred_ll_het',
                             'gt_phred_ll_homalt'])

        masks = ['False' if m is None or m.strip('(').strip(')') == 'False'
                 else m for m in self.family_masks]
        masks = [compiler.compile(m, m, 'eval') for m in masks]

        for gene, li in self.candidates():

            for row in li:
                cols = dict((col, row[col]) for col in req_cols)
                fams = [f for i, f in enumerate(self.families)
                        if masks[i] != 'False' and eval(masks[i], cols)]
                if len(fams) < self.args.min_kindreds: continue

                # an ordered dict.
                pdict = row.print_fields

                for fam in fams:
                    # populate with the fields required by the tools.
                    pdict["family_id"] = fam.family_id
                    pdict["family_members"] = ",".join("%s" % m for m in fam.subjects)
                    pdict["family_genotypes"] = ",".join([eval(str(s), cols) for s in fam.gts])
                    pdict["samples"] = ",".join([x.name or x.sample_id for x in fam.subjects if x.affected])
                    pdict["family_count"] = len(fams)
                    yield pdict

    def run(self):
        for i, s in enumerate(self.report_candidates()):
            if i == 0:
                print "\t".join(s.keys())
            print "\t".join(map(str, s.values()))


class AutoDom(GeminiInheritanceModel):
    model = "auto_dom"

    def candidates(self):
        for g, li in self.gen_candidates('gene'):
            yield g, li


class AutoRec(AutoDom):
    model = "auto_rec"


class DeNovo(GeminiInheritanceModel):
    model = "de_novo"

    def candidates(self):
        kins = self.args.min_kindreds
        for g, li in self.gen_candidates('gene' if kins is not None else None):
            yield g, li


class MendelViolations(GeminiInheritanceModel):
    model = "mendel_violations"

    def candidates(self):
        for g, li in self.gen_candidates(None):
            yield g, li
