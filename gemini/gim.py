#!/usr/bin/env python
import collections
import sys
from copy import copy
import re
import GeminiQuery
import sql_utils
import compiler
from gemini_constants import *
from .gemini_bcolz import filter, NoGTIndexException
from .mendelianerror import mendelian_error
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
        try:
            for i, family_id in enumerate(self.family_ids):
                gt_filter = self.family_masks[i]
                if isinstance(gt_filter, dict):
                    for k, flt in gt_filter.items():
                        variant_ids.update(filter(self.args.db, flt, {}))
                else:
                    if gt_filter == 'False': continue
                    # TODO: maybe we should just or these together and call filter once.
                    variant_ids.update(filter(self.args.db, gt_filter, {}))

            return sorted(set(variant_ids))
        except NoGTIndexException:
            return None

    def gen_candidates(self, group_key):
        if isinstance(group_key, basestring):
            group_key = op.itemgetter(group_key)

        q = self.query
        vids = self.bcolz_candidates()
        if vids is None:
            self.gq.run(q, needs_genotypes=True)

        elif len(vids) > 0:
            q = GeminiQuery.add_variant_ids_to_query(q, vids)
            self.gq.run(q, needs_genotypes=True)
        else:
            # no variants met the criteria
            raise StopIteration

        for grp_key, grp in it.groupby(self.gq, group_key):
            yield grp_key, grp

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
        kwargs = {'only_affected': not getattr(self.args, "allow_unaffected", False),
                  'strict': not self.args.lenient}
        if self.model == "mendel_violations":
            kwargs = {'only_affected': self.args.only_affected}
        for family in families:
            # e.g. family.auto_rec(gt_ll, min_depth)
            family_filter = getattr(family,
                    self.model)(gt_ll=self.args.gt_phred_ll,
                                min_depth=self.args.min_sample_depth,
                                **kwargs)

            self.family_masks.append(family_filter)
            self.family_ids.append(family.family_id)

    def report_candidates(self):
        args = self.args
        req_cols = ['gt_types', 'gts']
        if args.min_sample_depth and self.args.min_sample_depth > 0:
            req_cols.append('gt_depths')
        if args.gt_phred_ll or self.model == "mendel_violations":

            for col in ['gt_phred_ll_homref', 'gt_phred_ll_het',
                        'gt_phred_ll_homalt']:
                if col in self.gt_cols:
                    req_cols.append(col)

        is_mendel = False

        if isinstance(self.family_masks[0], dict):
            assert self.model == "mendel_violations"
            is_mendel = True
            masks = []
            # mdict contains filter for 'de novo', 'LOH', etc.
            for mdict in self.family_masks:
                m = {}
                for k, mask in mdict.items():
                    m[k] = 'False' if mask is None or mask.strip("(").strip(")") == 'False' else mask
                    if m[k] != 'False':
                        m[k] = compiler.compile(m[k], m[k], 'eval')

                masks.append(m)
        else:
            # 1 mask per family
            masks = ['False' if m is None or m.strip('(').strip(')') in
                     ('empty', 'False') else m for m in self.family_masks]
            masks = [compiler.compile(m, m, 'eval') if m != 'False' else 'False' for m in masks]

        requested_fams = None if not args.families else set(args.families.split(","))

        for gene, li in self.candidates():

            kindreds = set()
            to_yield = []

            seen = set()
            for row in li:
                # comp_het sets a family_id that met the filter. so we use it
                # for this check instead of checking all families.
                cur_fam = row.print_fields.get('family_id')

                cols = dict((col, row[col]) for col in req_cols)
                fams_to_test = enumerate(self.families) if cur_fam is None \
                                else [(i, f) for i, f in enumerate(self.families) if f.family_id == cur_fam]
                # limit to families requested by the user.
                if requested_fams is not None:
                    fams_to_test = ((i, f) for i, f in fams_to_test if f.family_id
                            in requested_fams)

                fams, models = [], []
                if is_mendel:
                    for i, f in fams_to_test:
                        mask_dict = masks[i]
                        for inh_model, mask in mask_dict.items():
                            if masks[i] != 'False' and eval(mask, cols):
                                if f in fams:
                                    models[-1] += ";" + inh_model
                                else:
                                    fams.append(f)
                                    models.append(inh_model)
                else:
                    fams = [f for i, f in fams_to_test
                            if masks[i] != 'False' and eval(masks[i], cols)]
                kindreds.update(f.family_id for f in fams)

                for fam in fams:
                    pdict = row.print_fields.copy()
                    kindreds.add(fam.family_id)
                    # get a *shallow* copy of the ordered dict.
                    # populate with the fields required by the tools.
                    pdict["family_id"] = fam.family_id
                    pdict["family_members"] = ",".join("%s" % m for m in fam.subjects)
                    pdict["family_genotypes"] = ",".join(eval(str(s), cols) for s in fam.gts)
                    pdict["samples"] = ",".join(x.name or x.sample_id for x in fam.subjects if x.affected)
                    pdict["family_count"] = len(fams)
                    if is_mendel:
                        pdict["violation"] = ";".join(models)
                        # TODO: check args, may need fam.subjects_with_parent
                        pdict['samples'] = fam.affecteds_with_parent if args.only_affected else fam.samples_with_parent
                        vs = []

                        if all(c in self.gt_cols for c in ('gt_phred_ll_homref', 'gt_phred_ll_het', 'gt_phred_ll_homalt')):
                            for s in pdict['samples']:
                                # mom, dad, kid
                                mdk = str(s.mom.genotype_lls + s.dad.genotype_lls + s.genotype_lls)
                                # get all 3 at once so we just do 1 eval.
                                vals = eval(mdk, cols)
                                vs.append(mendelian_error(vals[:3], vals[3:6], vals[6:], pls=True))

                            pdict["violation_prob"] = ",".join("%.5f" % v for v in vs)
                        pdict['samples'] = ",".join(s.name or str(s.sample_id) for s in pdict['samples'])

                    s = str(pdict)
                    if s in seen:
                        continue
                    seen.add(s)

                    to_yield.append(pdict)

            if len(kindreds) >= self.args.min_kindreds:
                for item in to_yield:
                    item['family_count'] = len(kindreds)
                    yield item

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


class CompoundHet(GeminiInheritanceModel):
    model = "comp_het"

    @property
    def query(self):
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

        if args.filter: query += " AND " + args.filter
        # we need to order results by gene so that we can sweep through the results
        return query + " ORDER BY gene"

    def _add_necessary_columns(self, custom_columns):
        """
        Convenience function to tack on columns that are necessary for
        the functionality of the tool but yet have not been specifically
        requested by the user.
        """
        # we need to add the variant's chrom, start and gene if
        # not already there.
        self.added = []
        for col in ("gene", "start", "alt", "variant_id"):
            if custom_columns.find(col) < 0:
                custom_columns += "," + col
                if col != "variant_id":
                    self.added.append(col)
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
                          comp_het_counter=[0]):
        """
        Refine candidate heterozygote pairs based on user's filters.
        """
        args = self.args
        # once we are in here, we know that we have a single gene.
        subjects_dict = self.subjects_dict

        candidates = samples_w_hetpair

        # TODO: we are filtering requested fams here before doing min-kindreds
        # count later. is this as expected?
        requested_fams = None if args.families is None else set(args.families.split(","))
        for idx, comp_het in enumerate(candidates):
            comp_het_counter[0] += 1
            #seen = set()
            for s in samples_w_hetpair[comp_het]:
                family_id = subjects_dict[s].family_id
                # NOTE: added this so each family only reported 1x.
                #if family_id in seen:
                #    continue
                if requested_fams is not None and not family_id in requested_fams:
                    continue
                #seen.add(family_id)

                ch_id = str(comp_het_counter[0])
                cid = "%s_%d_%d" % (ch_id, comp_het[0].row['variant_id'],
                                    comp_het[1].row['variant_id'])
                for i in (0, 1):
                    pdict = comp_het[i].row.print_fields.copy()
                    # set these to keep order in the ordered dict.
                    pdict["family_id"] = family_id
                    pdict["family_members"] = None
                    pdict["family_genotypes"] = None
                    pdict["samples"] = None
                    pdict["family_count"] = None
                    pdict["comp_het_id"] = cid

                    comp_het[i].row.print_fields = pdict
                    # TODO: check this yield.
                    yield comp_het[i].row

    def candidates(self):
        args = self.args
        idx_to_sample = self.gq.idx_to_sample

        from .family import Family
        self.gq._connect_to_database()
        fams = self.fams = Family.from_cursor(self.gq.c)

        self.subjects_dict = {}
        for f in fams:
            for s in fams[f].subjects:
                self.subjects_dict[s.name] = s

        for grp, li in self.gen_candidates('gene'):
            sample_hets = collections.defaultdict(lambda: collections.defaultdict(list))

            for row in li:

                gt_types, gt_bases, gt_phases = row['gt_types'], row['gts'], row['gt_phases']
                for famid, one_fam in fams.items():
                    # can phase each sample separately as they know their
                    # index into each gt array.
                    gt_phases, gt_bases = one_fam.famphase(gt_types, gt_phases, gt_bases, length_check=False)

                site = Site(row)
                # track each sample that is heteroyzgous at this site.
                for idx, gt_type in enumerate(gt_types):
                    if gt_type != HET:
                        continue
                    sample = idx_to_sample[idx]
                    # need to keep unaffecteds so we no if someone has the same
                    # pair.
                    #if args.only_affected and not self.subjects_dict[sample].affected:
                    #    continue
                    sample_site = copy(site)
                    sample_site.phased = gt_phases[idx]

                    if not sample_site.phased and not args.ignore_phasing:
                        continue

                    sample_site.gt = gt_bases[idx]
                    # add the site to the list of candidates for this sample/gene
                    sample_hets[sample][site.row['gene']].append(sample_site)

            # process the last gene seen
            samples_w_hetpair = self.find_valid_het_pairs(sample_hets)
            if not args.allow_unaffected:
                sd = self.subjects_dict
                # key is (site1, site2), value is list of samples
                # if a single unaffected sample shares this het pair, we remove
                # it from consideration.
                samples_w_hetpair = dict((site, samples) for site, samples in
                        samples_w_hetpair.items() if all(sd[s].affected for s in samples))

            yield grp, self.filter_candidates(samples_w_hetpair)

class Site(object):
    __slots__ = ('row', 'phased', 'gt')
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
