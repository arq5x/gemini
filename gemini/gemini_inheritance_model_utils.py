#!/usr/bin/env python
from collections import defaultdict
import GeminiQuery
import sql_utils
from mendelianerror import mendelian_error
from gemini_constants import *
import gemini_subjects as subjects
import itertools as it
import re

def get_prob(family_gt_lls, row):
    e = {}
    for k in ('gt_phred_ll_homref', 'gt_phred_ll_het', 'gt_phred_ll_homalt'):
        e[k] = row[k]

    for k, li in family_gt_lls.iteritems():
        e[k] = [eval(item, e) for item in li]
        # order is father, mother, child
    father = [e[k][0] for k in ("homref", "het", "homalt")]
    mother = [e[k][1] for k in ("homref", "het", "homalt")]
    child = [e[k][2] for k in ("homref", "het", "homalt")]
    return mendelian_error(father, mother, child, pls=True)

class GeminiInheritanceModelFactory(object):

    gt_cols = ('gts', 'gt_types', 'gt_phases', 'gt_depths', 'gt_ref_depths',
               'gt_alt_depths', 'gt_quals', 'gt_phred_ll_homref',
               'gt_phred_ll_het', 'gt_phred_ll_homalt')

    # https://github.com/arq5x/gemini/pull/436
    required_columns = ("family_id", "family_members",
                        "family_genotypes", "samples", "family_count")

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

    def report_candidates(self, candidates,
                          is_comp_het=False):
        """
        Print variants that meet the user's requirements
        If input is a tuple,
        """
        min_kindreds = self.args.min_kindreds
        family_count = len(candidates)

        if min_kindreds is not None and family_count < min_kindreds:
            return False

        candidate_keys = sorted(candidates.keys())
        fam_counts_by_variant = defaultdict(int)
        for (gene, family_id), li in candidates.items():
            for row in li:
                variant_id = row[3]
                fam_counts_by_variant[variant_id] += 1

        for (gene, family_id) in candidate_keys:
            for tup in candidates[(gene, family_id)]:

                # (row, family_gt_label, family_gt_cols) \
                violation = ""
                comp_het = ""
                if self.model == "mendel_violations":
                    (row, family_gt_label, family_gt_cols, variant_id,
                            family_gt_lls, violation) = tuple(tup)
                elif is_comp_het:
                    (row, family_gt_label, family_gt_cols, variant_id, comp_het,
                            subject) = tuple(tup)
                else:
                    (row, family_gt_label, family_gt_cols, family_id) = tuple(tup)

                e = {}
                for k in ('gt_types', 'gts', 'gt_depths'):
                    e[k] = row[k]

                if self.model == "mendel_violations":
                    prob = get_prob(family_gt_lls, row)
                    violation += ("\t%.3f" % prob)

                gt_cols = [str(eval(s, e)) for s in family_gt_cols]
                gt_lbls = family_gt_label

                if is_comp_het:
                    affected_samples = [subject]
                    for c in self.added:
                        if c in row.row:
                            row.row.pop(c)
                else:
                    affected_samples = self.get_samples(gt_cols, gt_lbls)

                print str(row) + "\t%s\t%s\t%s\t%s\t%i%s" % (
                             family_id,
                             ",".join(gt_lbls),
                             ",".join(gt_cols),
                             ",".join(affected_samples),
                             fam_counts_by_variant[variant_id],
                             ("\t" + violation + comp_het).rstrip())

    def get_samples(self, gts, lbls, splitter=re.compile("\/|\|"), model=None):
        """Hack..."""
        model = model or self.model
        if model == "auto_rec":
            # first will always be parent who must have ref
            ref = splitter.split(gts[0])[0]
            return [l.split("(", 1)[0] for gt, l in it.izip(gts, lbls)
                    if ref not in splitter.split(gt) and ";affected" in l]
        elif model == "auto_dom":
            # has to be the affected samples that are hets.
            return [l.split("(", 1)[0] for gt, l in it.izip(gts, lbls) if
                    ";affected" in l and len(set(splitter.split(gt))) > 1]
        elif model == "de_novo":
            # samples that have alleles that the parents do not.
            parents = set(splitter.split(gts[0]) + splitter.split(gts[1]))
            return [l.split("(", 1)[0] for gt, l in it.izip(gts[2:], lbls[2:]) if
                    set(splitter.split(gt)) - parents]
        elif model == "mendel_violations":
            mom, dad = set(splitter.split(gts[0])), set(splitter.split(gts[1]))
            # denovo
            samps = self.get_samples(gts, lbls, model="de_novo")
            # u disomy. kid didn't get one from each parent. and LOH
            samps += [l.split("(", 1)[0] for gt, l in it.izip(gts[2:], lbls[2:])
                      if not (set(splitter.split(gt)) - mom
                              and set(splitter.split(gt)) - dad)]
            return samps

    def _cull_families(self):
        """
        If the user has asked to restric the analysis to a specific set
        of families, then we need to prune the list of possible families
        to that specific subset.
        """
        pass

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
        self.family_gt_phred_lls = []
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
                if self.model == "mendel_violations":
                    self.family_gt_phred_lls.append(family.get_genotype_lls())

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
        if self.model in ("auto_rec", "auto_dom") or \
           (self.model == "de_novo" and self.args.min_kindreds is not None):

            # we require the "gene" column for the auto_* tools
            self.query = sql_utils.ensure_columns(self.query, ['gene'])
            if self.args.filter:
                self.query += " AND gene is not NULL ORDER BY chrom, gene"
            else:
                self.query += " WHERE gene is not NULL ORDER BY chrom, gene"
        self.query = sql_utils.ensure_columns(self.query, ['variant_id'])

    def get_header(self, gqh, is_comp_het=False):
        h = "\t".join(self.required_columns)

        header = gqh + "\t" + h
        if self.model == "mendel_violations":
            return header + "\tviolation\tviolation_prob"
        elif is_comp_het:
            gqh = "\t".join([hd.strip() for hd in gqh.split("\t") if not hd.strip() in self.added])
            return gqh + "\t" + h + "\tcomp_het_id"
        else:
            return header

    def _get_gene_only_candidates(self):
        """
        Identify candidates that meet the user's criteria AND affect genes.
        """
        # collect family info
        self._get_family_info()

        # run the query applying any genotype filters provided by the user.
        self._construct_query()
        self.gq.run(self.query, needs_genotypes=True)

        print self.get_header(self.gq.header)

        # yield the resulting variants for this familiy
        self.candidates = defaultdict(list)
        prev_gene = None
        for row in self.gq:

            curr_gene = row['gene']
            variant_id = row.row['variant_id']

            # report any candidates for the previous gene
            if curr_gene != prev_gene and prev_gene is not None:
                self.report_candidates(self.candidates)
                # reset for the next gene
                self.candidates = defaultdict(list)

            # test the variant for each family in the db
            for idx, family_id in enumerate(self.family_ids):
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
                if self.model == "mendel_violations":
                    family_gt_phred_lls = self.family_gt_phred_lls[idx]
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
                key = (curr_gene, family_id)
                self.candidates[key].append([row,
                                             family_gt_labels,
                                             family_gt_cols,
                                             variant_id])
                if self.model == "mendel_violations":
                    self.candidates[key][-1].append(family_gt_phred_lls)
                    self.candidates[key][-1].append(",".join(violations))

            prev_gene = curr_gene

        # report any candidates for the last gene
        self.report_candidates(self.candidates)

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

        print self.get_header(self.gq.header)
        gene = None

        for row in self.gq:
            #variant_id = row.row.pop('variant_id')
            variant_id = row.row['variant_id']

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
                if self.model == "mendel_violations":
                    family_gt_phred_lls = self.family_gt_phred_lls[idx]
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
                key = (gene, family_id)
                candidates[key].append([row, family_gt_labels, family_gt_cols, variant_id])

                if self.model == "mendel_violations":
                    candidates[key][-1].append(family_gt_phred_lls)
                    candidates[key][-1].append(",".join(violations))
            self.report_candidates(candidates)

