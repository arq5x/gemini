#!/usr/bin/env python
from collections import defaultdict
import GeminiQuery
import sql_utils
from gemini_constants import *
import gemini_subjects as subjects


class GeminiInheritanceModelFactory(object):

    gt_cols = ('gts', 'gt_types', 'gt_phases', 'gt_depths', 'gt_ref_depths',
               'gt_alt_depths', 'gt_quals', 'gt_phred_ll_homref',
               'gt_phred_ll_het', 'gt_phred_ll_homalt')

    # https://github.com/arq5x/gemini/pull/436
    required_columns = ("variant_id", "family_id", "family_members",
                        "family_genotypes", "affected_samples", "family_count")

    def __init__(self, args, model):

        # default to all genotype columns and all columns
        if not args.columns:
            args.columns = "*," + ", ".join(self.gt_cols)

        self.args = args
        self.model = model
        self.gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)

    def get_candidates(self):
        """
        Report candidate variants that meet the requested inheritance model.
        """
        if self.model in ["auto_dom", "auto_rec"] \
           or (self.model == "de_novo" and self.args.min_kindreds is not None):
            self._get_gene_only_candidates()
        else:
            self._get_all_candidates()

    @classmethod
    def report_candidates(self, candidates, is_violation=False,
                          min_kindreds=1):
        """
        Print variants that meet the user's requirements
        If input is a tuple,
        """
        family_count = len(candidates)

        if family_count < min_kindreds:
            return False

        candidate_keys = sorted(candidates.keys())
        for (variant_id, gene) in candidate_keys:
            for tup in candidates[(variant_id, gene)]:

                # (row, family_gt_label, family_gt_cols) \
                if is_violation:
                    (row, family_gt_label, family_gt_cols, family_id, violation) = tuple(tup)
                    violation += "\t"
                else:
                    (row, family_gt_label, family_gt_cols, family_id) = tuple(tup)
                    violation = ""

                e = {}
                for k in ('gt_types', 'gts', 'gt_depths'):
                    e[k] = row[k]

                affected_samples = [x.split("(")[0] for x in family_gt_label if ";affected" in x]
                print str(row) + "\tvariant_id:%s\t%s%s\t%s\t%s\t%s\t%i" % (variant_id, violation, family_id,
                                              ",".join(str(s) for s in family_gt_label),
                                              ",".join(str(eval(s, e)) for s
                                                  in family_gt_cols),
                                              ",".join(affected_samples),
                                              family_count)

    def _cull_families(self):
        """
        If the user has asked to restric the analysis to a specific set
        of families, then we need to prune the list of possible families
        to that specific subset.
        """

    def _get_family_info(self):
        """
        Extract the relevant genotype filters, as well all labels
        for each family in the database.
        """
        families = subjects.get_families(self.args.db, self.args.families)
        self.family_ids = []
        self.family_masks = []
        self.family_gt_labels = []
        self.family_gt_columns = []
        self.family_dp_columns = []
        for family in families:

            family_filter = None

            if self.model == "auto_rec":
                family_filter = family.get_auto_recessive_filter(gt_ll=self.args.gt_phred_ll)
            elif self.model == "auto_dom":
                family_filter = family.get_auto_dominant_filter(gt_ll=self.args.gt_phred_ll)
            elif self.model == "de_novo":
                family_filter = family.get_de_novo_filter(self.args.only_affected, gt_ll=self.args.gt_phred_ll)
            elif self.model == "mendel_violations":
                family_filter = family.get_mendelian_violation_filter(gt_ll=self.args.gt_phred_ll)

            if family_filter != "False" and family_filter is not None:
                self.family_masks.append(family_filter)
                self.family_gt_labels.append(family.get_genotype_labels())
                self.family_gt_columns.append(family.get_genotype_columns())
                self.family_dp_columns.append(family.get_genotype_depths())
                self.family_ids.append(family.family_id)

    def _construct_query(self):
        """
        Construct the relevant query based on the user's requests.
        """
        if self.args.columns is not None:
            # the user only wants to report a subset of the columns
            self.query = "SELECT " + str(self.args.columns) + " FROM variants "
        else:
            # report the kitchen sink
            self.query = "SELECT chrom, start, end, * \
                    , gts, gt_types, gt_phases, gt_depths, \
                    gt_ref_depths, gt_alt_depths, gt_quals, \
                    gt_phred_ll_homref, gt_phred_ll_het, gt_phred_ll_homalt \
                    FROM variants "

        # add any non-genotype column limits to the where clause
        if self.args.filter:
            self.query += " WHERE " + self.args.filter

        # auto_rec and auto_dom candidates should be limited to
        # variants affecting genes.
        if self.model == "auto_rec" or self.model == "auto_dom"\
        or (self.model == "de_novo" and self.args.min_kindreds is not None):

            # we require the "gene" column for the auto_* tools
            self.query = sql_utils.ensure_columns(self.query, ['gene'])
            if self.args.filter:
                self.query += " AND gene is not NULL ORDER BY chrom, gene"
            else:
                self.query += " WHERE gene is not NULL ORDER BY chrom, gene"
        self.query = sql_utils.ensure_columns(self.query, ['variant_id'])

    def get_header(self, is_violation_query):
        h = "\t".join(self.required_columns)
        gqh = self.gq.header
        # strip variant_id as they didn't request it, but we added it for the
        # required columns
        if gqh.endswith("\tvariant_id"):
            gqh, _ = gqh.rsplit("\t", 1)
        if is_violation_query:
            return gqh + "\tviolation\t" + h
        else:
            return gqh + "\t" + h

    def _get_gene_only_candidates(self):
        """
        Identify candidates that meet the user's criteria AND affect genes.
        """
        # collect family info
        self._get_family_info()

        # run the query applying any genotype filters provided by the user.
        self._construct_query()
        self.gq.run(self.query, needs_genotypes=True)

        is_violation_query = isinstance(self.family_masks[0], dict)
        print self.get_header(is_violation_query)

        # yield the resulting variants for this familiy
        self.candidates = defaultdict(list)
        prev_gene = None
        for row in self.gq:

            curr_gene = row['gene']
            variant_id = row.row.pop('variant_id')

            # report any candidates for the previous gene
            if curr_gene != prev_gene and prev_gene is not None:
                self.report_candidates(self.candidates, is_violation_query,
                        self.args.min_kindreds)
                # reset for the next gene
                self.candidates = defaultdict(list)

            # test the variant for each family in the db
            for idx, fam_id in enumerate(self.family_ids):
                family_genotype_mask = self.family_masks[idx]
                family_gt_labels = self.family_gt_labels[idx]
                family_gt_cols = self.family_gt_columns[idx]
                family_dp_cols = self.family_dp_columns[idx]

                # interrogate the genotypes present in each family member to
                # conforming to the genetic model being tested
                e = {}
                for c in ('gt_types', 'gts', 'gt_depths', 'gt_phred_ll_homalt',
                          'gt_phred_ll_het', 'gt_phred_ll_homref'):
                    e[c] = row[c]

                # skip if the variant doesn't meet a recessive model
                # for this family
                violations = []
                if is_violation_query:
                    for violation, mask in family_genotype_mask.items():
                        if eval(mask, e):
                            violations.append(violation)
                    if len(violations) == 0:
                        continue
                else:
                    if not eval(family_genotype_mask, e):
                        continue

                # make sure each sample's genotype had sufficient coverage.
                # otherwise, ignore
                insufficient_depth = False
                for col in family_dp_cols:
                    depth = int(eval(col, e))
                    if depth < self.args.min_sample_depth:
                        insufficient_depth = True
                        break
                if insufficient_depth is True:
                    continue

                # if it meets a recessive model, add it to the list
                # of candidates for this gene.
                self.candidates[(variant_id, curr_gene)].append([row,
                                                        family_gt_labels,
                                                        family_gt_cols,
                                                        fam_id])
                if is_violation_query:
                    self.candidates[(variant_id, curr_gene)][-1].append(",".join(violations))

            prev_gene = curr_gene

        # report any candidates for the last gene
        self.report_candidates(self.candidates, is_violation_query,
                self.args.min_kindreds)

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
        self.gq.run(self.query, needs_genotypes=True)

        is_violation_query = isinstance(self.family_masks[0], dict)
        print self.get_header(is_violation_query)
        gene = None

        for row in self.gq:
            variant_id = row.row.pop('variant_id')
            candidates = defaultdict(list)
            cols = {}
            for col in self.gt_cols:
                cols[col] = row[col]

            # test the variant for each family in the db
            for idx, family_id in enumerate(self.family_ids):
                family_genotype_mask = self.family_masks[idx]
                family_gt_labels = self.family_gt_labels[idx]
                family_gt_cols = self.family_gt_columns[idx]
                family_dp_cols = self.family_dp_columns[idx]

                # interrogate the genotypes present in each family member to
                # conforming to the genetic model being tested

                # skip if the variant doesn't meet a recessive model
                # for this family
                violations = []
                if is_violation_query:
                    for violation, mask in family_genotype_mask.items():
                        if eval(mask, cols):
                            violations.append(violation)
                    if len(violations) == 0:
                        continue
                else:
                    if not eval(family_genotype_mask, cols):
                        continue

                # make sure each sample's genotype had sufficient coverage.
                # otherwise, ignore
                insufficient_depth = False
                for col in family_dp_cols:
                    depth = int(eval(col, cols))
                    if depth < self.args.min_sample_depth:
                        insufficient_depth = True
                        break
                if insufficient_depth is True:
                    continue

                # shoe-horn the variant so we can use report_candidates.
                key = (variant_id, gene)
                candidates[key].append([row, family_gt_labels, family_gt_cols, family_id])
                if is_violation_query:
                    candidates[key].append(",".join(violations))

            self.report_candidates(candidates, is_violation_query,
                    self.args.min_kindreds)
