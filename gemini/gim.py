#!/usr/bin/env python
import sys
from collections import Counter, defaultdict
import GeminiQuery
import sql_utils
import compiler
from gemini_constants import *
from .gemini_bcolz import filter, NoGTIndexException
from .mendelianerror import mendelian_error
import itertools as it
import operator as op
from inheritance import Family


class GeminiInheritanceModel(object):

    required_columns = ("gene", "family_id", "family_members",
                        "family_genotypes", "samples", "family_count")

    def __init__(self, args):

        self.args = args
        self.gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
        self.added = []

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

        query = sql_utils.ensure_columns(query, ['variant_id', 'gene'])
        # add any non-genotype column limits to the where clause
        if self.args.filter:
            query += " WHERE " + self.args.filter

        if hasattr(self.args, 'X'):
            if self.args.X == []:
                self.args.X = ['chrX', 'X']
            part = "chrom IN (%s)" % ", ".join("'%s'" % x for x in self.args.X)
            if " WHERE " in query:
                query += " AND " + part
            else:
                query += " WHERE " + part

        # auto_rec and auto_dom candidates should be limited to
        # variants affecting genes.
        if self.model in ("auto_rec", "auto_dom", "comp_het") or \
           (self.model == "de_novo" and self.args.min_kindreds is not None):

            # we require the "gene" column for the auto_* tools
            if " WHERE " in query:
                query += " AND gene is not NULL"
            else:
                query += " WHERE gene is not NULL"

        # only need to order by gene for comp_het and when min_kindreds is used.
        if self.model == "comp_het" or not (
                self.args.min_kindreds in (None, 1)
                and (self.args.families is None
                     or not "," in self.args.families)):
            query += " ORDER by chrom, gene"

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
        self.families = families = Family.from_cursor(self.gq.conn).values()
        args = self.args

        self.family_ids = []
        self.family_masks = []
        kwargs = {'only_affected': not getattr(self.args, "allow_unaffected", False),
                  'min_gq': args.min_gq}
        if self.model == "mendel_violations":
            kwargs = {'only_affected': self.args.only_affected}
        if self.model != "comp_het" and self.model != "mendel_violations":
            if hasattr(self.args, 'lenient'):
                kwargs['strict'] = not self.args.lenient
        elif self.model == "comp_het":
            kwargs['pattern_only'] = self.args.pattern_only
        if hasattr(self.args, 'gt_phred_ll'):
            kwargs['gt_ll'] = self.args.gt_phred_ll

        if self.model in ('x_rec', 'x_dom', 'x_denovo'):
            kwargs.pop('only_affected')

        requested_fams = None if not args.families else set(args.families.split(","))

        for family in families:
            if requested_fams is None or family.family_id in requested_fams:
                # e.g. family.auto_rec(gt_ll, min_depth)
                family_filter = getattr(family, self.model)(
                                    min_depth=self.args.min_sample_depth,
                                    **kwargs)
            else:
                family_filter = 'False'

            self.family_masks.append(family_filter)
            self.family_ids.append(family.family_id)

    def report_candidates(self):
        args = self.args
        req_cols = ['gt_types', 'gts']
        if args.min_sample_depth and self.args.min_sample_depth > 0:
            req_cols.append('gt_depths')
        if getattr(args, 'gt_phred_ll', False) or self.model == "mendel_violations":

            for col in ['gt_phred_ll_homref', 'gt_phred_ll_het',
                        'gt_phred_ll_homalt']:
                if col in self.gt_cols:
                    req_cols.append(col)
        if args.min_gq != 0 and 'gt_quals' in self.gt_cols:
            req_cols.append('gt_quals')

        is_mendel = False

        if any(isinstance(m, dict) for m in self.family_masks):
            assert self.model == "mendel_violations"
            is_mendel = True
            masks = []
            # mdict contains filter for 'de novo', 'LOH', etc.
            for mdict in self.family_masks:
                if isinstance(mdict, basestring):
                    masks.append(mdict)
                    continue
                m = {}
                for k, mask in mdict.items():
                    m[k] = 'False' if mask is None or mask.strip("(").strip(")") == 'False' else mask
                    if m[k] != 'False':
                        m[k] = compiler.compile(m[k], m[k], 'eval')

                masks.append(m)
        else:
            # 1 mask per family
            masks = ['False' if m is None or m is False or m.strip('(').strip(')') in
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
                for col in self.added:
                    try:
                        del row.print_fields[col]
                    except ValueError:
                        pass

                cols = dict((col, row[col]) for col in req_cols)
                fams_to_test = enumerate(self.families) if cur_fam is None \
                                else [(i, f) for i, f in enumerate(self.families) if f.family_id == cur_fam]
                # limit to families requested by the user.
                if requested_fams is not None:
                    fams_to_test = ((i, f) for i, f in fams_to_test if f.family_id in requested_fams)

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
                            pdict["violation_prob"] = ""
                            for s in pdict['samples']:
                                # mom, dad, kid
                                mdk = str(s.mom.genotype_lls + s.dad.genotype_lls + s.genotype_lls)
                                # get all 3 at once so we just do 1 eval.
                                vals = eval(mdk, cols)
                                vs.append(mendelian_error(vals[:3], vals[3:6], vals[6:], pls=True))

                            pdict["violation_prob"] = ",".join("%.5f" % v for v in vs if v is not None)
                        pdict['samples'] = ",".join(s.name or str(s.sample_id) for s in pdict['samples'])

                    s = str(pdict)
                    if s in seen:
                        continue
                    seen.add(s)

                    to_yield.append(pdict)

            if 0 < len(kindreds) >= self.args.min_kindreds:
                if 'comp_het_id' in to_yield[0]:
                    counts = Counter((item['comp_het_id'], item['family_id']) for item in to_yield)
                    # remove singletons.
                    to_yield = [item for item in to_yield if counts[(item['comp_het_id'], item['family_id'])] > 1]
                    # re-check min_kindreds
                    if len(set(item['family_id'] for item in to_yield)) < self.args.min_kindreds:
                        continue

                for item in to_yield:
                    item['family_count'] = len(kindreds)
                    yield item

    def run(self):
        has_gts = False
        from .gemini_bcolz import gt_cols_types
        for i, s in enumerate(self.report_candidates()):
            if i == 0:
                has_gts = [x[0] for x in gt_cols_types if x[0] in s] or False
                print "\t".join(s.keys())
            if has_gts:
                for col in has_gts:
                    s[col] = str(s[col]).replace('\n', '')
            print "\t".join(map(str, s.values()))


class XRec(GeminiInheritanceModel):
    model = "x_rec"

    def candidates(self):
        for g, li in self.gen_candidates('gene'):
            yield g, li

class XDenovo(XRec):
    model = "x_denovo"

class XDom(XRec):
    model = "x_dom"


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
        return query + " ORDER BY chrom, gene"

    def _add_necessary_columns(self, custom_columns):
        """
        Convenience function to tack on columns that are necessary for
        the functionality of the tool but yet have not been specifically
        requested by the user.
        """
        # we need to add the variant's chrom, start and gene if
        # not already there.
        custom_columns = [x.strip() for x in custom_columns.split(",")]
        if "*" in custom_columns:
            return ",".join(custom_columns)
        self.added = []
        for col in ("gene", "chrom", "start", "ref", "alt", "variant_id"):
            if not col in custom_columns:
                custom_columns.append(col)
                if col != "variant_id":
                    self.added.append(col)
        return ",".join(custom_columns)

    def filter_candidates(self, candidates,
                          comp_het_counter=[0]):
        """
        Refine candidate heterozygote pairs based on user's filters.
        """
        args = self.args
        # once we are in here, we know that we have a single gene.

        requested_fams = None if args.families is None else set(args.families.split(","))
        for idx, comp_het in enumerate(candidates):
            comp_het_counter[0] += 1

            for fam_ch in candidates[comp_het]:
                if fam_ch['priority'] > args.max_priority:
                    continue
                # when to use affected_unphased?
                for subject in (fam_ch['candidates'] if args.pattern_only else fam_ch['affected_phased'] + fam_ch['affected_unphased']):
                    family_id = subject.family_id

                    if requested_fams is not None and not family_id in requested_fams:
                        continue

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
                        pdict['priority'] = fam_ch['priority']

                        comp_het[i].row.print_fields = pdict
                        # TODO: check this yield.
                        yield comp_het[i].row

    def candidates(self):
        args = self.args

        self.gq._connect_to_database()
        fams = self.fams = Family.from_cursor(self.gq.conn)

        if args.families:
            fams = {f: fam for f, fam in fams.items() if f in set(args.families.split(","))}

        for grp, li in self.gen_candidates('gene'):
            samples_w_hetpair = defaultdict(list)
            sites = []
            for row in li:

                gt_types, gt_bases, gt_phases = row['gt_types'], row['gts'], row['gt_phases']
                site = Site(row)
                site.gt_phases, site.gt_bases, site.gt_types = gt_phases, gt_bases, gt_types
                sites.append(site)

            for i, site1 in enumerate(sites[:-1], start=1):
                for site2 in sites[i:]:

                    for family_id, fam in fams.items():

                        ch = fam.comp_het_pair(site1.gt_types, site1.gt_bases,
                                               site2.gt_types, site2.gt_bases,
                                               site1.gt_phases, site2.gt_phases,
                                               ref1=site1.row['ref'],
                                               alt1=site1.row['alt'],
                                               ref2=site2.row['ref'],
                                               alt2=site2.row['alt'],
                                               allow_unaffected=args.allow_unaffected,
                                               fast_mode=True,
                                               pattern_only=args.pattern_only)

                        if not ch['candidate']: continue

                        samples_w_hetpair[(site1, site2)].append(ch)
            yield grp, self.filter_candidates(samples_w_hetpair)

class Site(object):
    __slots__ = ('row', 'gt_phases', 'gt_bases', 'gt_types')
    def __init__(self, row):
        self.row = row
        self.gt_phases = None
        self.gt_bases = None
        self.gt_types = None

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
