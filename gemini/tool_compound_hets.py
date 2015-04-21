#!/usr/bin/env python
import os
import sys
import collections
import re
from copy import copy

import GeminiQuery
import gemini_utils as util
from gemini_constants import *
from gemini_inheritance_model_utils import GeminiInheritanceModelFactory as Factory
import gemini_subjects as subjects
from itertools import groupby
from operator import itemgetter

class Site(object):
    def __init__(self, row):
        self.row = row
        self.phased = None
        self.gt = None

    def __eq__(self, other):
        return self.row['chrom'] == other.row['chrom'] and \
               self.row['start'] == other.row['start']

    def __repr__(self):
        return ",".join([self.row['chrom'],
                         str(self.row['start']),
                         str(self.row['end'])])

    def __hash__(self):
        "hash the site based on chrom+start"
        return sum(ord(c) for c in self.row['chrom']) + int(self.row['start'])




class CompoundHet(Factory):


    def create_query(self):
        """
        Construct a query to identify candidate compound heterozygotes
        based on the user's columns and filters
        """
        args = self.args
        if args.columns is not None:
            custom_columns = self._add_necessary_columns(str(args.columns))
            query = "SELECT " + custom_columns + \
                    " FROM variants " + \
                    " WHERE (is_exonic = 1 or impact_severity != 'LOW') "
        else:
            # report the kitchen sink
            query = "SELECT *" + \
                    ", gts, gt_types, gt_phases, gt_depths, \
                    gt_ref_depths, gt_alt_depths, gt_quals" + \
                    " FROM variants " + \
                    " WHERE (is_exonic = 1 or impact_severity != 'LOW') "

        # add any non-genotype column limits to the where clause
        if args.filter:
            query += " AND " + args.filter

        # we need to order results by gene so that we can sweep through the results
        query += " ORDER BY gene"
        return query

    def _add_necessary_columns(self, custom_columns):
        """
        Convenience function to tack on columns that are necessary for
        the functionality of the tool but yet have not been specifically
        requested by the user.
        """
        # we need to add the variant's chrom, start and gene if
        # not already there.
        self.added = []
        if custom_columns.find("gene") < 0:
            custom_columns += ", gene"
            self.added.append("gene")
        if custom_columns.find("start") < 0:
            custom_columns += ", start"
            self.added.append("start")
        if custom_columns.find("alt") < 0:
            custom_columns += ", alt"
            self.added.append("alt")
        if custom_columns.find("variant_id") < 0:
            custom_columns += ", variant_id"

        return custom_columns

    def find_valid_het_pairs(self, sample_hets):
        """
        Identify candidate heterozygote pairs.
        """
        args = self.args
        samples_w_hetpair = collections.defaultdict(list)
        splitter = re.compile("\||/")
        for sample in sample_hets:
            for gene in sample_hets[sample]:

                # we only care about combinations, not permutations
                # (e.g. only need site1,site2, not site1,site2 _and site2,site1)
                # thus we can do this in a ~ linear pass instead of a ~ N^2 pass
                for idx, site1 in enumerate(sample_hets[sample][gene]):
                    for site2 in sample_hets[sample][gene][idx + 1:]:

                        # expand the genotypes for this sample at each site into
                        # it's composite alleles.  e.g. A|G -> ['A', 'G']
                        alleles_site1 = []
                        alleles_site2 = []
                        if not args.ignore_phasing:
                            alleles_site1 = site1.gt.split('|')
                            alleles_site2 = site2.gt.split('|')
                        else:
                            # split on phased (|) or unphased (/) genotypes
                            alleles_site1 = splitter.split(site1.gt)
                            alleles_site2 = splitter.split(site2.gt)

                        # it is only a true compound heterozygote IFF
                        # the alternates are on opposite haplotypes.
                        if not args.ignore_phasing:
                            # return the haplotype on which the alternate allele
                            # was observed for this sample at each candidate het.
                            # site. e.g., if ALT=G and alleles_site1=['A', 'G']
                            # then alt_hap_1 = 1.  if ALT=A, then alt_hap_1 = 0
                            if "," in str(site1.row['alt']) or \
                               "," in str(site2.row['alt']):
                                sys.stderr.write("WARNING: Skipping candidate for sample"
                                                 " %s b/c variants with mult. alt."
                                                 " alleles are not yet supported. The sites are:"
                                                 " %s and %s.\n" % (sample, site1, site2))
                                continue

                            alt_hap_1 = alleles_site1.index(site1.row['alt'])
                            alt_hap_2 = alleles_site2.index(site2.row['alt'])

                        # Keep as a candidate if
                        #   1. phasing is considered AND the alt alleles are on
                        #      different haplotypes
                        #   2. the user doesn't care about phasing.
                        # TODO: Phase based on parental genotypes.
                        if (not args.ignore_phasing and alt_hap_1 != alt_hap_2) \
                            or args.ignore_phasing:
                            samples_w_hetpair[(site1,site2)].append(sample)

        return samples_w_hetpair


    def filter_candidates(self, samples_w_hetpair,
                          family_gt_labels,
                          family_gt_cols,
                          comp_het_counter=[0]):
        """
        Refine candidate heterozygote pairs based on user's filters.
        """
        args = self.args
        # eliminate comp_hets with unaffected individuals if
        # only affected individuals are required.
        # once we are in here, we know that we have a single gene.

        candidates = {}
        if args.only_affected:
            for comp_het in samples_w_hetpair:
                num_affected = sum(self.subjects_dict[s].affected \
                                    for s in samples_w_hetpair[comp_het])
                if num_affected == len(samples_w_hetpair[comp_het]):
                    candidates[comp_het] = samples_w_hetpair[comp_het]
        else:
            candidates = samples_w_hetpair

        # catalog the set of families that have a comp_het in this gene
        family_count = collections.Counter()
        for comp_het in candidates:
            for s in samples_w_hetpair[comp_het]:
                family_id = self.subjects_dict[s].family_id
                family_count[family_id] += 1

        # were there enough families with a compound het in this gene?
        # keys of (variant_id, gene) vals of [row, family_gt_label, family_gt_cols,
        # family_id, comp_het_id]
        filtered_candidates = collections.defaultdict(list)
        if len(family_count) >= args.min_kindreds:
            for idx, comp_het in enumerate(candidates):
                comp_het_counter[0] += 1
                for s in samples_w_hetpair[comp_het]:
                    family_id = self.subjects_dict[s].family_id
                    if args.families is not None and family_id not in args.families.split(','):
                        continue

                    ch_id = str(comp_het_counter[0])
                    for i in (0, 1):
                        row = comp_het[i].row

                        filtered_candidates[(row['gene'], family_id)]\
                            .append((row,
                                     family_gt_labels[family_id],
                                     family_gt_cols[family_id],
                                     row['variant_id'],
                                     ch_id,
                                     s))
            self.report_candidates(filtered_candidates, is_comp_het=True)

    def get_candidates(self):
        self.get_compound_hets()

    def get_compound_hets(self):
        """
        Report candidate compound heterozygotes.
        """
        args = self.args
        gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
        idx_to_sample = gq.idx_to_sample
        self.subjects_dict = subjects.get_subjects(args)

        # run the query applying any genotype filters provided by the user.
        gq.run(self.create_query())

        families = subjects.get_families(args.db, args.families)
        family_gt_labels, family_gt_cols = {}, {}
        for family in families:
            family_gt_labels[family.family_id] = family.get_genotype_labels()
            family_gt_cols[family.family_id] = family.get_genotype_columns()

        # output header
        print self.get_header(gq.header, is_comp_het=True)

        # Collect all of the genic heterozygotes for each sample / gene
        for gene, row_list in groupby(gq, itemgetter('gene')):
            sample_hets = collections.defaultdict(lambda: collections.defaultdict(list))
            for row in row_list:

                gt_types, gt_bases, gt_phases = row['gt_types'], row['gts'], row['gt_phases']
                site = Site(row)
                # track each sample that is heteroyzgous at this site.
                for idx, gt_type in enumerate(gt_types):
                    if gt_type != HET:
                        continue
                    sample = idx_to_sample[idx]
                    sample_site = copy(site)
                    sample_site.phased = gt_phases[idx]

                    if not sample_site.phased and not args.ignore_phasing:
                        continue

                    sample_site.gt = gt_bases[idx]
                    # add the site to the list of candidates for this sample/gene
                    sample_hets[sample][site.row['gene']].append(sample_site)

            # process the last gene seen
            samples_w_hetpair = self.find_valid_het_pairs(sample_hets)
            self.filter_candidates(samples_w_hetpair,
                              family_gt_labels,
                              family_gt_cols)

def run(parser, args):
    if os.path.exists(args.db):
        CompoundHet(args, "comp_het").get_compound_hets()

