#!/usr/bin/env python
import sqlite3
import os
import sys
from collections import defaultdict
import GeminiQuery
from gemini_constants import *
import gemini_subjects as subjects


class GeminiInheritanceModelFactory(object):

    def __init__(self, args, model):
        self.args = args
        self.model = model
        self.gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)


    def get_candidates(self):
        """
        Report candidate variants that meet the requested inheritance model.
        """
        if self.model in ["auto_dom", "auto_rec"]:
            self._get_gene_only_candidates()
        else:
            self._get_all_candidates()


    def _report_candidates(self):
        """
        Print variants that meet the user's requirements
        """
        num_families = self.candidates.keys()

        if len(num_families) >= self.args.min_kindreds:

            for (gene, family_id) in self.candidates:
                for (row, family_gt_label, family_gt_cols) \
                    in self.candidates[(gene, family_id)]:

                    gt_types = row['gt_types']
                    gts = row['gts']

                    print str(family_id) + "\t" + \
                        ",".join([str(s) for s in family_gt_label]) + \
                        "\t", \
                        ",".join([str(eval(s)) for s in family_gt_cols]) + \
                        "\t",
                    print row

    def _get_family_info(self):
        """
        Extract the relevant genotype filters, as well all labels
        for each family in the database.
        """
        families = subjects.get_families(self.gq.c)
        self.family_ids = []
        self.family_masks = []
        self.family_gt_labels = []
        self.family_gt_columns = []
        for family in families:

            family_filter = None
            
            if self.model == "auto_rec":
                family_filter = family.get_auto_recessive_filter()
            elif self.model == "auto_dom":
                family_filter = family.get_auto_dominant_filter()
            elif self.model == "de_novo":
                family_filter = family.get_de_novo_filter()

            if family_filter != "False" and family_filter is not None:
                self.family_masks.append(family_filter)
                self.family_gt_labels.append(family.get_genotype_labels())
                self.family_gt_columns.append(family.get_genotype_columns())
                self.family_ids.append(family.family_id)

    def _construct_query(self):
        """
        Construct the relevant query based on the user's requests.
        """
        if self.args.columns is not None:
            # the user only wants to report a subset of the columns
            self.query = "SELECT " + str(self.args.columns) + \
                    " FROM variants "
        else:
            # report the kitchen sink
            self.query = "SELECT *" + \
                    ", gts, gt_types, gt_phases, gt_depths, \
                    gt_ref_depths, gt_alt_depths, gt_quals" + \
                    " FROM variants "

        # add any non-genotype column limits to the where clause
        if self.args.filter:
            self.query += " WHERE " + self.args.filter
        
        # auto_rec and auto_dom candidates should be limited to
        # variants affecting genes.
        if self.model == "auto_rec" or self.model == "auto_dom":
            if self.args.filter:
                self.query += " AND gene is not NULL ORDER BY gene"
            else:
                self.query += " WHERE gene is not NULL ORDER BY gene"

    def _get_gene_only_candidates(self):
        """
        Identify candidates that meet the user's criteria AND affect genes.
        """
        # collect family info
        self._get_family_info()

        # run the query applying any genotype filters provided by the user.
        self._construct_query()
        self.gq.run(self.query)

        # print a header
        print "family_id\tfamily_members\tfamily_genotypes\t",
        print self.gq.header

        # yield the resulting variants for this familiy
        self.candidates = defaultdict(list)
        prev_gene = None
        for row in self.gq:

            curr_gene = row['gene']
        
            # report any candidates for the previous gene
            if curr_gene != prev_gene and prev_gene is not None:
                self._report_candidates()
                # reset for the next gene
                self.candidates = defaultdict(list)

            # test the variant for each family in the db
            for idx, fam_id in enumerate(self.family_ids):
                family_genotype_mask = self.family_masks[idx]
                family_gt_labels = self.family_gt_labels[idx]
                family_gt_cols = self.family_gt_columns[idx]

                # interrogate the genotypes present in each family member to 
                # conforming to the genetic model being tested
                gt_types = row['gt_types']
                gts = row['gts']

                # skip if the variant doesn't meet a recessive model
                # for this family
                if not eval(family_genotype_mask):
                    continue
                    

                # if it meets a recessive model, add it to the list
                # of candidates for this gene.
                self.candidates[(curr_gene, fam_id)].append((row,
                                                        family_gt_labels, 
                                                        family_gt_cols))

            prev_gene = curr_gene
    
        # report any candidates for the last gene
        self._report_candidates()
        
    def _get_all_candidates(self):
        """
        Identify candidates that meet the user's criteria no matter where
        they occur in the genome.
        """
        """
        Identify candidates that meet the user's criteria AND affect genes.
        """
        # collect family info
        self._get_family_info()

        # run the query applying any genotype filters provided by the user.
        self._construct_query()
        self.gq.run(self.query)

        # print a header
        print "family_id\tfamily_members\tfamily_genotypes\t",
        print self.gq.header

        for row in self.gq:

            # test the variant for each family in the db
            for idx, fam_id in enumerate(self.family_ids):
                family_genotype_mask = self.family_masks[idx]
                family_gt_labels = self.family_gt_labels[idx]
                family_gt_cols = self.family_gt_columns[idx]

                # interrogate the genotypes present in each family member to 
                # conforming to the genetic model being tested
                gt_types = row['gt_types']
                gts = row['gts']

                # skip if the variant doesn't meet a recessive model
                # for this family
                if not eval(family_genotype_mask):
                    continue

                print str(family_id) + "\t" + \
                    ",".join([str(s) for s in family_gt_label]) + \
                    "\t", \
                    ",".join([str(eval(s)) for s in family_gt_cols]) + \
                    "\t",
                print row