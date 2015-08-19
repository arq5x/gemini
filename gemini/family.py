"""
Create filters for given inheritance models.
See: https://github.com/arq5x/gemini/issues/388
"""
from __future__ import print_function
import sys

from collections import defaultdict, Counter
import itertools as it
import operator as op
import re
try:
    reduce
except NameError:
    from functools import reduce

HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)

valid_gts = (
    'gts',
    'gt_types',
    'gt_phases',
    'gt_depths',
    'gt_ref_depths',
    'gt_alt_depths',
    'gt_quals',
    'gt_copy_numbers',
    'gt_phred_ll_homref',
    'gt_phred_ll_het',
    'gt_phred_ll_homalt',
)

def _bracket(other):
    o = str(other)
    if o in ("True", "False"): return o
    return "(%s)" % o

class ostr(str):
    def __and__(self, other):
        if other is None: return self
        if other is empty: return self
        if self is empty: return other
        return ostr("%s and %s" % (_bracket(self), _bracket(other)))

    def __or__(self, other):
        if other is None: return self
        if other is empty: return self
        if self is empty: return other
        return ostr("%s or %s" % (_bracket(self), _bracket(other)))

    def __nonzero__(self):
        raise Exception("shouldn't be here. use & instead of 'and'. and wrap in parens")

empty = ostr("empty")

class Sample(object):
    """
    >>> a, b = Sample(1, True), Sample(2, True)
    >>> a.gt_types == HOM_REF
    'gt_types[0] == HOM_REF'

    >>> a.gt_types == b.gt_types
    'gt_types[0] == gt_types[1]'

    >>> a.gt_types != b.gt_types
    'gt_types[0] != gt_types[1]'

    >>> (a.gt_phred_ll_homref < 2) & (a.gt_phred_ll_homalt > 2)
    '(gt_phred_ll_homref[0] < 2) and (gt_phred_ll_homalt[0] > 2)'

    >>> (a.gt_phred_ll_homref < 2) | (b.gt_phred_ll_homalt > 2)
    '(gt_phred_ll_homref[0] < 2) or (gt_phred_ll_homalt[1] > 2)'
    """

    __slots__ = ('sample_id', 'name', 'affected', 'sex', 'gender', 'mom', 'dad',
                 'family_id', '_i')

    def __init__(self, sample_id, affected, gender=None, name=None,
                 family_id=None):
        #assert isinstance(sample_id, (long, int)), sample_id
        assert affected in (True, False, None)
        self.sample_id = sample_id
        self.name = name or sample_id
        self.affected = affected
        self.mom = None
        self.dad = None
        self.sex = self.gender = gender
        self.family_id = family_id
        # _i is used to maintain the order in which they came in.
        self._i = None

    @property
    def genotype_lls(self):
        return [self.gt_phred_ll_homref,
                self.gt_phred_ll_het,
                self.gt_phred_ll_homalt]

    def __getattr__(self, gt_field):
        assert gt_field in valid_gts, gt_field
        return Filter(self.sample_id, gt_field)

    def __repr__(self):
        c = self.__class__.__name__
        s = "%s(%s" % (c, self.name or self.sample_id)
        s += (";affected" if self.affected else (";unaffected"
                  if self.affected is False else ";unknown"))
        if self.gender is not None:
            s += ";%s" % self.gender
        return s + ")"

    def __str__(self):
        r = repr(self).split("(", 1)[1]
        return "%s(%s" % (self.name, r)


class Filter(object):
    def __init__(self, sample_id, gt_field):
        self.gt_field = gt_field
        if isinstance(sample_id, (long, int)):
            self.sample0 = sample_id - 1
        else:
            self.sample0 = sample_id

    def in_(self, li):
        return reduce(op.or_, [self == i for i in li])

    def __lt__(self, o):
        return ostr("%s[%s] < %s" % (self.gt_field, self.sample0, o))

    def __le__(self, o):
        return ostr("%s[%s] <= %s" % (self.gt_field, self.sample0, o))

    def __gt__(self, o):
        return ostr("%s[%s] > %s" % (self.gt_field, self.sample0, o))

    def __ge__(self, o):
        return ostr("%s[%s] >= %s" % (self.gt_field, self.sample0, o))

    def __eq__(self, o):
        return ostr("%s[%s] == %s" % (self.gt_field, self.sample0, o))

    def __ne__(self, o):
        return ostr("%s[%s] != %s" % (self.gt_field, self.sample0, o))

    def __str__(self):
        return ostr("%s[%s]" % (self.gt_field, self.sample0))

    def __repr__(self):
        return ostr("%s[%s]" % (self.gt_field, self.sample0))

    def __and__(self, other):
        raise Exception("shouldn't be here. wrap &/| statements in parens")

    __or__ = __and__

    def __nonzero__(self):
        raise Exception("shouldn't be here. use & instead of 'and'")

class Family(object):
    """
    >>> mom = Sample(1, affected=False)
    >>> dad = Sample(2, affected=False)
    >>> kid = Sample(3, affected=True)
    >>> Family([mom, dad, kid], 'a').auto_rec()
    '(gt_types[2] == HOM_ALT) and ((gt_types[0] != HOM_ALT) and (gt_types[1] != HOM_ALT))'

    # unknowns don't count
    >>> kid2 = Sample(4, affected=None)

    >>> Family([mom, dad, kid, kid2], 'a').auto_rec()
    '(gt_types[2] == HOM_ALT) and ((gt_types[0] != HOM_ALT) and (gt_types[1] != HOM_ALT))'

    >>> Family([mom, dad, kid, kid2], 'a').auto_rec(gt_ll=1, min_depth=10)
    '(((gt_types[2] == HOM_ALT) and (gt_phred_ll_homalt[2] <= 1)) and (((gt_types[0] != HOM_ALT) and (gt_types[1] != HOM_ALT)) and ((gt_phred_ll_homalt[0] > 1) and (gt_phred_ll_homalt[1] > 1)))) and ((gt_depths[2] >= 10) and ((gt_depths[0] >= 10) and (gt_depths[1] >= 10)))'

    """

    def __init__(self, subjects, fam_id):
        assert len(set(s.sample_id for s in subjects)) == len(subjects), subjects
        self.subjects = subjects
        self.family_id = fam_id
        # mostly for testing. this gets set when reading from db.
        if all(s._i is None for s in self.subjects):
            for i, s in enumerate(self.subjects):
                s._i = i

    @classmethod
    def from_ped(klass, ped):
        """
        return a dict keyed by family_id with parent/kid relations defined from a ped file.
        """
        def agen():
            for toks in (l.rstrip().split() for l in open(ped) if l[0] != "#"):
                toks.append(toks[1]) # name
                yield toks
        return klass._from_gen(agen())

    def __len__(self):
        return len(self.subjects)

    def famphase(self, gt_types, gt_phases, gt_bases,
                 length_check=True,
                 _splitter=re.compile("\||/")):
        """
        >>> HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)
        >>> f = Family([Sample("mom", False), Sample("dad", False),
        ...             Sample("kid", True)], "fam")
        >>> f.subjects[2].mom = f.subjects[0]
        >>> f.subjects[2].dad = f.subjects[1]
        >>> f.famphase([HOM_REF, HOM_ALT, HET],
        ...            [False, False, False],
        ...            ["A/A", "A/C", "C/A"])
        ([False, False, True], ['A/A', 'A/C', 'A|C'])

        >>> f.famphase([HOM_REF, HOM_REF, HET],
        ...            [False, False, False],
        ...            ["A/A", "A/C", "C/A"])
        ([False, False, False], ['A/A', 'A/C', 'C/A'])

        >>> f.famphase([HOM_REF, HET, HET],
        ...            [False, False, False],
        ...            ["A/A", "A/C", "C/A"])
        ([False, False, True], ['A/A', 'A/C', 'A|C'])

        >>> f.famphase([HET, HET, HET],
        ...            [False, False, False],
        ...            ["A/C", "A/C", "A/C"])
        ([False, False, False], ['A/C', 'A/C', 'A/C'])

        >>> f.famphase([HOM_REF, HET, HET],
        ...            [False, False, False],
        ...            ["AA/A", "AA/C", "C/AA"])
        ([False, False, True],  ['AA/A', 'AA/C', 'AA|C'])

        >>> f.famphase([HOM_REF, HET, HET],
        ...            [False, False, False],
        ...            ['G/G', 'AA/C', 'A/C'])
        ([False, False, False],  ['G/G', 'AA/C', 'A/C'])

        >>> f.famphase([HOM_REF, HET, HET],
        ...            [False, False, False],
        ...            ['G/G', 'A/C', 'A/C'])
        ([False, False, False],  ['G/G', 'A/C', 'A/C'])
        """

        # NOTE: this modifies in-place
        # subjects are in same order as gt_types and _i is the index.
        HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)
        if length_check:
            assert len(self.subjects) == len(gt_types) == len(gt_bases)
        assert isinstance(gt_bases[0], basestring)
        assert isinstance(gt_types[0], int) or int(gt_types[0]) == gt_types[0]

        for s in (subj for subj in self.subjects if gt_types[subj._i] == HET):
            # here we have a phaseable (HET) person ..
            if s.mom is None or s.dad is None: continue
            # ... with a mom and dad
            ## cant be same.
            if gt_types[s.mom._i] == gt_types[s.dad._i]: continue
            ## cant have unknown
            if UNKNOWN in (gt_types[s.mom._i], gt_types[s.dad._i]): continue

            kid_bases = set(_splitter.split(gt_bases[s._i]))
            mom_bases = _splitter.split(gt_bases[s.mom._i])
            dad_bases = _splitter.split(gt_bases[s.dad._i])

            parent_bases = set(mom_bases + dad_bases)
            # can't phase kid with de-novo

            if kid_bases - parent_bases:
                sys.stderr.write("skipping variant due to apparent de_novo in kid\n")
                continue

            # no alleles from dad
            if len(kid_bases - set(dad_bases)) == len(kid_bases):
                sys.stderr.write("skipping variant due to no alleles from dad (apparent mendelian error)\n")
                continue

            if len(kid_bases - set(mom_bases)) == len(kid_bases):
                sys.stderr.write("skipping variant due to no alleles from mom (apparent mendelian error)\n")
                continue

            # should be able to phase here
            if gt_types[s.mom._i] in (HOM_REF, HOM_ALT):
                assert gt_types[s.mom._i] in (HOM_REF, HOM_ALT)
                mom_allele = mom_bases[0]
                dad_alleles = dad_bases
                dad_allele = next(d for d in dad_alleles if d != mom_allele)


                gt_bases[s._i] = "%s|%s" % (mom_allele, dad_allele)
            else:
                assert gt_types[s.dad._i] in (HOM_REF, HOM_ALT)

                dad_allele = dad_bases[0]
                mom_alleles = mom_bases
                mom_allele = next(m for m in mom_alleles if m != dad_allele)

                gt_bases[s._i] = "%s|%s" % (mom_allele, dad_allele)

            gt_phases[s._i] = True

        return gt_phases, gt_bases

    def to_ped(self, fh=sys.stdout, header=True):
        if header:
            fh.write("#family_id sample_id paternal_id maternal_id sex phenotype\n")
        for s in self.subjects:
            paternal_id = (s.dad and s.dad.name) or "-9"
            maternal_id = (s.mom and s.mom.name) or "-9"
            phenotype = {True: '2', False: '1'}.get(s.affected, "-9")
            fh.write(" ".join((str(self.family_id), s.name, paternal_id, maternal_id, s.sex or '-9', phenotype)))
            fh.write("\n")

    @classmethod
    def from_cursor(klass, cursor):
        keys = "sample_id|family_id|name|paternal_id|maternal_id|sex|phenotype".split("|")
        def agen():
            for row in cursor.execute("select %s from samples" % ",".join(keys)):
                if not isinstance(row, dict):
                    row = dict(zip(keys, row))
                yield (row['family_id'], row['sample_id'], row['paternal_id'],
                       row['maternal_id'], row['sex'], str(row['phenotype']),
                       row['name'])
        return klass._from_gen(agen())

    @classmethod
    def _from_gen(klass, gen):
        fams = defaultdict(dict)
        pheno_lookup = {'1': False, '2': True}
        gender_lookup = {'1': 'male', '2': 'female'}
        for i, (fam_id, indv, pat_id, mat_id, sex, pheno, name) in enumerate(gen):

            assert indv not in fams[fam_id]
            s = fams[fam_id][name] = Sample(indv, pheno_lookup.get(pheno),
                                            gender=gender_lookup.get(sex))
            s.mom = mat_id
            s.dad = pat_id
            # name in gemini is actually the id from the ped.
            # sample_id in gemini si the primary key
            s.name = name
            s.family_id = fam_id
            s._i = i

        ofams = {}
        for fam_id, fam_dict in fams.items():
            ofams[fam_id] = []
            for name in fam_dict:
                sample = fam_dict[name]
                # convert sample_id to dad or None
                sample.dad = fam_dict.get(sample.dad)
                sample.mom = fam_dict.get(sample.mom)
                ofams[fam_id].append(sample)

            ofams[fam_id] = Family(ofams[fam_id], fam_id)
            # maintain the order in which they came in.
            ofams[fam_id].subjects.sort(key=op.attrgetter('_i'))
        return ofams

    def __repr__(self):
        return "%s([%s])" % (self.__class__.__name__,
                           ", ".join(repr(s) for s in self.subjects))

    @property
    def affecteds(self):
        return [s for s in self.subjects if s.affected]

    @property
    def affecteds_with_parent(self):
        return [s for s in self.affecteds if not None in (s.mom, s.dad)]

    @property
    def samples_with_parent(self):
        return [s for s in self.subjects if not None in (s.mom, s.dad)]

    @property
    def unaffecteds(self):
        return [s for s in self.subjects if s.affected is False]

    @property
    def unknown(self):
        return [s for s in self.subjects if s.affected is None]

    # so we can do, e.g. fam.gts and get the list of required strings.
    def __getattr__(self, gt_field):
        assert gt_field in valid_gts, gt_field
        return [getattr(s, gt_field) for s in self.subjects]

    def _restrict_to_min_depth(self, min_depth, unknowns=False):
        if min_depth is not None and min_depth > 0:
            af = reduce(op.and_, [s.gt_depths >= min_depth for s in self.affecteds], empty)
            un = reduce(op.and_, [s.gt_depths >= min_depth for s in self.unaffecteds], empty)
            if unknowns:
                return af & un & reduce(op.and_, [s.gt_depths >= min_depth for s in self.unknown], empty)
            return af & un
        else:
            return None

    def auto_dom(self, min_depth=0, gt_ll=False, strict=True, only_affected=True):
        """
        If strict then all affected kids must have at least 1 affected parent.
        parent.
        Parents of affected can't have unknown phenotype (for at least 1 kid)
        """
        if len(self.affecteds) == 0:
            sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
            if strict:
                return 'False'
        af = reduce(op.and_, [s.gt_types == HET for s in self.affecteds], empty)
        if len(self.unaffecteds) and only_affected:
            un = reduce(op.and_, [(s.gt_types != HET) & (s.gt_types != HOM_ALT) for s in self.unaffecteds])
        else:
            un = None
        depth = self._restrict_to_min_depth(min_depth)
        if gt_ll:
            af &= reduce(op.and_, [s.gt_phred_ll_het <= gt_ll for s in self.affecteds], empty)
            if len(self.unaffecteds) and only_affected:
                un &= reduce(op.and_, [s.gt_phred_ll_het > gt_ll for s in self.unaffecteds])
        # need at least 1 kid with parent who has the mutation
        # parents can't have unkown phenotype.
        kid_with_known_parents = False
        # all affected kids must have at least 1 affected parent (or no parents)
        kid_with_parents = False
        for kid in self.affecteds:
            # if they have no parents, don't require it
            if kid.mom is None and kid.dad is None:
                continue
            # mom or dad must be affected.
            kid_with_parents = True
            if strict and not any(p is not None and p.affected for p in (kid.mom, kid.dad)):
                return 'False'
            # parents can't have unknown phenotype.
            if (kid.mom and kid.dad):
                if (kid.mom.affected is not None) and (kid.dad.affected is not None):
                    kid_with_known_parents = True
                # if he has a mom and dad that arent unknown, at least one of them must be affected
                if not None in (kid.mom.affected, kid.dad.affected):
                    if not (kid.mom.affected or kid.dad.affected): return 'False'

        if strict and not kid_with_known_parents:
            return 'False'

        if not kid_with_parents:
            if len(self.affecteds) > 0:
                sys.stderr.write("WARNING: using affected without parents for family %s for autosomal dominant test. Use strict to prevent this.\n" % self.family_id)
        return af & un & depth

    def auto_rec(self, min_depth=0, gt_ll=False, strict=True, only_affected=True):
        """
        If strict, then if parents exist, they must be het for all affecteds
        """
        if len(self.affecteds) == 0:
            sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
            return 'False'
        af = reduce(op.and_, [s.gt_types == HOM_ALT for s in self.affecteds])
        if only_affected:
            un = reduce(op.and_, [s.gt_types != HOM_ALT for s in
                                  self.unaffecteds], empty)
        else:
            un = None
        if strict:
            # if parents exist, they must be het or affected for all affecteds
            # if both parents are not het then it's a de novo.
            usable_kid = False
            for kid in self.affecteds:
                usable_kid = usable_kid or (kid.mom and kid.dad)
                for parent in (kid.mom, kid.dad):
                    if parent is not None:
                        af &= parent.gt_types == HET
                        if parent.affected:
                            sys.stderr.write("WARNING: auto-recessive called on family "
                                    "%s where affected has affected parents\n" % self.family_id)
                            return "False"
                if not usable_kid:
                    sys.stderr.write("WARNING: auto-recessive called on family "
                            "%s where no affected has parents\n" % self.family_id)

        depth = self._restrict_to_min_depth(min_depth)
        if gt_ll:
            af &= reduce(op.and_, [s.gt_phred_ll_homalt <= gt_ll for s in self.affecteds])
            if only_affected:
                un &= reduce(op.and_, [s.gt_phred_ll_homalt > gt_ll for s in self.unaffecteds])

        return af & un & depth

    def de_novo(self, min_depth=0, gt_ll=False, strict=True, only_affected=True):
        """
        all affected must be het.
        all unaffected must be homref.
        if strict, all affected kids must have unaffected parents.
        >>> f = Family([Sample("mom", False), Sample("dad", False),
        ...             Sample("kid", True)], "fam")
        >>> f.gt_types = [HOM_ALT, HOM_ALT, HET]
        >>> f.subjects[2].mom = f.subjects[0]
        >>> f.subjects[2].dad = f.subjects[1]

        #>>> f.de_novo()
        #'False'

        >>> f2 = Family([Sample("mom", False), Sample("dad", False),
        ...              Sample("kid", True), Sample("kid2", False)], "fam")
        >>> f2.subjects[2].mom = f.subjects[0]
        >>> f2.subjects[2].dad = f.subjects[1]
        >>> f2.subjects[3].mom = f.subjects[0]
        >>> f2.subjects[3].dad = f.subjects[1]

        >>> f2.gt_types = [HOM_REF, HOM_REF, HET, HET]

        #>>> f2.de_novo()
        #'False'

        #>>> f2.de_novo(only_affected=False, strict=False)
        #'True'

        """
        if len(self.affecteds) == 0:
            sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
            return 'False'
        af = reduce(op.and_, [s.gt_types == HET for s in self.affecteds])
        un = empty
        if only_affected:
            un = reduce(op.and_, [s.gt_types == HOM_REF for s in self.unaffecteds], empty)
        if gt_ll:
            af &= reduce(op.and_, [s.gt_phred_ll_het <= gt_ll for s in self.affecteds], empty)
            if only_affected:
                un &= reduce(op.and_, [s.gt_phred_ll_homref <= gt_ll for s in self.unaffecteds], empty)

        if only_affected:
            un2 = reduce(op.and_, [s.gt_types == HOM_ALT for s in
                         self.unaffecteds], empty)
            if gt_ll:
                un2 &= reduce(op.and_, [s.gt_phred_ll_homalt <= gt_ll for s in self.unaffecteds], empty)
            un |= un2

        # at least 1 affected kid must have unaffected parents
        un_parents = False
        for kid in self.affecteds:
            if kid.mom and kid.mom.affected is False and kid.dad and kid.dad.affected is False:
                un_parents = True
            # can't have a parent with the variant
            for parent in (kid.mom, kid.dad):
                if parent is None: continue
                un &= (parent.gt_types == HOM_REF) | (parent.gt_types == HOM_ALT)
            if kid.mom and kid.dad:
                # otherwise, they could be HOM_REF, HOM_ALT and a het kid is not
                # de_novo
                un &= (kid.mom.gt_types == kid.dad.gt_types)

        if not un_parents:
            return "False"

        depth = self._restrict_to_min_depth(min_depth)
        if strict:
            # if a parent is affected it's not de novo.
            for kid in self.affecteds:
                for parent in (kid.mom, kid.dad):
                    if parent is not None and parent.affected:
                        return 'False'
        return af & un & depth

    def mendel_plausible_denovo(self, min_depth=0, gt_ll=False, only_affected=False):
        """
        kid == HET and dad, mom == HOM_REF or dad, mom == HOM_ALT
        only use kids with both parents present.
        """
        if only_affected:
            subset = self.affecteds_with_parent
        else:
            subset = self.samples_with_parent
        if len(subset) == 0: return 'False'

        depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
        # join exprs by or. so just look for any kid that meets these within a family.
        exprs = []
        for s in subset:
            # kid het, parents hom_ref
            expr = (s.gt_types == HET) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HOM_REF)
            if gt_ll:
                expr &= (s.gt_phred_ll_het <= gt_ll) & (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_homref <= gt_ll)
            exprs.append(expr)
            # kid het, parents hom_alt
            expr = (s.gt_types == HET) & (s.mom.gt_types == HOM_ALT) & (s.dad.gt_types == HOM_ALT)
            if gt_ll:
                expr &= (s.gt_phred_ll_het <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
            exprs.append(expr)
        return reduce(op.or_, exprs) & depth

    def mendel_implausible_denovo(self, min_depth=0, gt_ll=False, only_affected=False):
        # everyone is homozygote. kid is opposit of parents.
        # only use kids with both parents present.
        if only_affected:
            subset = self.affecteds_with_parent
        else:
            subset = self.samples_with_parent
        if len(subset) == 0: return 'False'

        depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
        exprs = []
        for s in subset:
            # kid hom_alt, parents hom_ref
            expr = (s.gt_types == HOM_ALT) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HOM_REF)
            if gt_ll:
                expr &= (s.gt_phred_ll_homalt <= gt_ll) & (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_homref <= gt_ll)
            exprs.append(expr)
            # parents hom_alt kid homref
            expr = (s.gt_types == HOM_REF) & (s.mom.gt_types == HOM_ALT) & (s.dad.gt_types == HOM_ALT)
            if gt_ll:
                expr &= (s.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
            exprs.append(expr)
        return reduce(op.or_, exprs) & depth

    def mendel_uniparental_disomy(self, min_depth=0, gt_ll=False,
                                  only_affected=False):
        # parents are opposite homs, kid matches one of them (but should be
        # het).
        if only_affected:
            subset = self.affecteds_with_parent
        else:
            subset = self.samples_with_parent
        if len(subset) == 0: return 'False'
        depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)

        # join exprs with or
        exprs = []
        for s in subset:
            for gtkid in (HOM_REF, HOM_ALT):
                # mom homref, dad hom_alt.
                expr = (s.gt_types == gtkid) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HOM_ALT)
                if gt_ll:
                    if gtkid == HOM_REF:
                        expr &= (s.gt_phred_ll_homref <= gt_ll)
                    else:
                        expr &= (s.gt_phred_ll_homalt <= gt_ll)
                    expr &= (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
                exprs.append(expr)

                # mom homalt, dad hom_ref.
                expr = (s.gt_types == gtkid) & (s.dad.gt_types == HOM_REF) & (s.mom.gt_types == HOM_ALT)
                if gt_ll:
                    if gtkid == HOM_REF:
                        expr &= (s.gt_phred_ll_homref <= gt_ll)
                    else:
                        expr &= (s.gt_phred_ll_homalt <= gt_ll)
                    expr &= (s.dad.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll)
                exprs.append(expr)

        return reduce(op.or_, exprs) & depth

    def mendel_LOH(self, min_depth=0, gt_ll=False, only_affected=False):
        # kid and one parent are opposite homozygotes other parent is het.
        if only_affected:
            subset = self.affecteds_with_parent
        else:
            subset = self.samples_with_parent
        if len(subset) == 0: return 'False'

        depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
        exprs = []  # joined with or
        for s in subset:
            # kid hom_alt, mom hom_ref, dad het.
            e = (s.gt_types == HOM_ALT) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HET)
            if gt_ll:
                e &= (s.gt_phred_ll_homalt <= gt_ll) & (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_het <= gt_ll)
            exprs.append(e)

            # kid hom_ref, mom het, dad hom_alt
            e = (s.gt_types == HOM_REF) & (s.mom.gt_types == HET) & (s.dad.gt_types == HOM_ALT)
            if gt_ll:
                e &= (s.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_het <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
            exprs.append(e)

            # kid hom_alt, mom het, dad hom_ref
            e = (s.gt_types == HOM_ALT) & (s.mom.gt_types == HET) & (s.dad.gt_types == HOM_REF)
            if gt_ll:
                e &= (s.gt_phred_ll_homalt <= gt_ll) & (s.mom.gt_phred_ll_het <= gt_ll) & (s.dad.gt_phred_ll_homref <= gt_ll)
            exprs.append(e)

            # kid hom_ref, mom hom_alt, dad het
            e = (s.gt_types == HOM_REF) & (s.mom.gt_types == HOM_ALT) & (s.dad.gt_types == HET)
            if gt_ll:
                e &= (s.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll) & (s.dad.gt_phred_ll_het <= gt_ll)
            exprs.append(e)

        return reduce(op.or_, exprs) & depth

    def mendel_violations(self, min_depth=0, gt_ll=False, only_affected=False):
        """
        >>> f = Family([Sample("mom", False), Sample("dad", False),
        ...             Sample("kid", True)], "fam")
        >>> f.subjects[2].mom = f.subjects[0]
        >>> f.subjects[2].dad = f.subjects[1]
        >>> r = f.mendel_violations()
        >>> for k in r:
        ...     print(k)
        ...     print(r[k])
        uniparental disomy
        (((((gt_types[kid] == HOM_REF) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_ALT)) or (((gt_types[kid] == HOM_REF) and (gt_types[dad] == HOM_REF)) and (gt_types[mom] == HOM_ALT))) or (((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_ALT))) or (((gt_types[kid] == HOM_ALT) and (gt_types[dad] == HOM_REF)) and (gt_types[mom] == HOM_ALT))
        plausible de novo
        (((gt_types[kid] == HET) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_REF)) or (((gt_types[kid] == HET) and (gt_types[mom] == HOM_ALT)) and (gt_types[dad] == HOM_ALT))
        loss of heterozygosity
        (((((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HET)) or (((gt_types[kid] == HOM_REF) and (gt_types[mom] == HET)) and (gt_types[dad] == HOM_ALT))) or (((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HET)) and (gt_types[dad] == HOM_REF))) or (((gt_types[kid] == HOM_REF) and (gt_types[mom] == HOM_ALT)) and (gt_types[dad] == HET))
        implausible de novo
        (((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_REF)) or (((gt_types[kid] == HOM_REF) and (gt_types[mom] == HOM_ALT)) and (gt_types[dad] == HOM_ALT))
        """
        return {'plausible de novo': self.mendel_plausible_denovo(min_depth,
                                                                  gt_ll,
                                                                  only_affected),
                'implausible de novo': self.mendel_implausible_denovo(min_depth,
                                                                      gt_ll,
                                                                      only_affected),
                'uniparental disomy': self.mendel_uniparental_disomy(min_depth,
                                                                     gt_ll,
                                                                     only_affected),
                'loss of heterozygosity': self.mendel_LOH(min_depth, gt_ll,
                                                          only_affected)
                }

    def _get_ref_alt(self, gt_types, gt_bases,
                     _splitter=re.compile("\||/")):
        """
        Guess the ref and alt. Mostly for convenience for comp_het functions,
        as we should know these anyway.
        """
        ref, alt = None, None
        for i, gt in enumerate(gt_types):
            if gt == HOM_REF:
                ref = _splitter.split(gt_bases[i])[0]
            elif gt == HOM_ALT:
                alt = _splitter.split(gt_bases[i])[0]
            elif "/" in gt_bases[i]:
                _ref, _alt = gt_bases[i].split("/")
                if ref is None:
                    ref = _ref
                if alt is None and _ref != _alt:
                    alt = _alt
        # fall back to allele frequency
        if ref is None or alt is None or ref == alt:
            c = Counter()
            for b in gt_bases:
                c.update(_splitter.split(b))
            if ref is None:
                ref = c.most_common(1)[0][0]
                if ref == alt:
                    ref = c.most_common(2)[1][0]
            if alt is None:
                alt = c.most_common(2)[1][0]
                if ref == alt:
                    alt = c.most_common(1)[0][0]
        return ref, alt

    def _comp_het_pair_pattern(self,
                               gt_types1, gt_nums1,
                               gt_types2, gt_nums2,
                               gt_phases1, gt_phases2):
        """
        + kid has to be phased het at both sites.
        + kid has to have alts on different chroms.
        + neither parent can be hom_alt at either site.
        + if either parent is phased at both sites and matches the kid, exclude.
        + if either parent is het at both sites, priority is reduced
        """

        # already phased before sending here.
        ret = {'candidates': [], 'priority': 4}
        for kid in self.samples_with_parent:
            if gt_nums1[kid._i] == gt_nums2[kid._i]: continue
            if not (gt_types1[kid._i] == HET and gt_types2[kid._i] == HET): continue
            #if not (gt_phases1[kid._i] and gt_phases2[kid._i]): continue
            if gt_types1[kid.mom._i] == HOM_ALT or gt_types2[kid.dad._i] == HOM_ALT: continue
            mom, dad = kid.mom, kid.dad

            kid_phased = gt_phases1[kid._i] and gt_phases2[kid._i]
            dad_phased = gt_phases1[dad._i] and gt_phases2[dad._i]
            mom_phased = gt_phases1[mom._i] and gt_phases2[mom._i]

            if kid_phased and dad_phased and (gt_nums1[dad._i] == gt_nums1[kid._i]) and (gt_nums2[dad._i] == gt_nums2[kid._i]):
                continue
            if kid_phased and mom_phased and (gt_nums1[mom._i] == gt_nums1[kid._i]) and (gt_nums2[mom._i] == gt_nums2[kid._i]):
                continue

            if kid_phased and dad_phased and mom_phased and gt_types1[dad._i] != gt_types2[dad._i] and gt_types1[mom._i] != gt_types2[mom._i]:
                priority = 1

            elif kid_phased and gt_types1[dad._i] != gt_types1[mom._i] and gt_types2[dad._i] != gt_types2[mom._i]:
                # parents are unphased hets at different sites.
                priority = 1
            else:
                priority = 2
                for parent in (kid.mom, kid.dad):
                    # unphased het
                    if gt_types2[parent._i] == gt_types1[parent._i] == HET:
                        priority += 1

            ret['candidates'].append(kid)
            ret['priority'] = min(ret['priority'], priority)
        ret['candidate'] = len(ret['candidates']) > 0
        return ret


    def comp_het_pair(self, gt_types1, gt_bases1,
                      gt_types2, gt_bases2,
                      gt_phases1=None,
                      gt_phases2=None,
                      ref1=None, alt1=None,
                      ref2=None, alt2=None,
                      allow_unaffected=False,
                      pattern_only=False,
                      _splitter=re.compile("\||/")):
        """
        Each of the sites here must have passed the comp_het() filter.
        This further checks that a give pair is comp_het.

        if pattern_only is False, affected/unaffected status is ignored.
        """
        if gt_phases1 is None:
            gt_phases1 = ["|" in b for b in gt_bases1]
        if gt_phases2 is None:
            gt_phases2 = ["|" in b for b in gt_bases2]

        if ref1 is None and alt1 is None:
            ref1, alt1 = self._get_ref_alt(gt_types1, gt_bases1)

        if ref2 is None and alt2 is None:
            ref2, alt2 = self._get_ref_alt(gt_types2, gt_bases2)

        self.famphase(gt_types1, gt_phases1, gt_bases1,
                      length_check=False)
        self.famphase(gt_types2, gt_phases2, gt_bases2,
                      length_check=False)

        gt_bases1 = [_splitter.split(b) for b in gt_bases1]
        gt_bases2 = [_splitter.split(b) for b in gt_bases2]

        # get in (0, 1) format instead of (A, T)
        ra = [ref1, alt1, "."]
        gt_nums1 = [(ra.index(b[0]), ra.index(b[1])) for b in gt_bases1]
        ra = [ref2, alt2, "."]
        gt_nums2 = [(ra.index(b[0]), ra.index(b[1])) for b in gt_bases2]

        if pattern_only:
            return self._comp_het_pair_pattern(gt_types1, gt_nums1,
                                               gt_types2, gt_nums2,
                                               gt_phases1, gt_phases2)

        for un in self.unaffecteds:
            if gt_types2[un._i] == HOM_ALT or gt_types1[un._i] == HOM_ALT:
                return {'candidate': False}

        ret = {'affected_phased': [], 'unaffected_phased': [],
               'unaffected_unphased': [], 'affected_unphased': [],
               'affected_skipped': [], 'candidates': []}

        aff = None
        for aff in self.affecteds:
            if gt_types1[aff._i] != HET or gt_types2[aff._i] != HET:
                ret['affected_skipped'].append(aff)
                # Remove candidates where an affected from the same family does
                # NOT share the same het pair.
                ret['candidate'] = False
                continue

            aff_phased = gt_phases1[aff._i] and gt_phases2[aff._i]
            # on same chrom.
            if aff_phased and gt_nums1[aff._i] == gt_nums2[aff._i]:
                ret['affected_skipped'].append(aff)
                # Remove candidates where an affected from the same family does
                # NOT share the same het pair.
                ret['candidate'] = False
                continue

            if not 'candidate' in ret: ret['candidate'] = True
            if aff_phased:
                ret['affected_phased'].append(aff)
            else:
                ret['affected_unphased'].append(aff)
            ret['candidates'].append(aff)

        del aff
        for un in self.unaffecteds:
            if gt_types1[un._i] != HET or gt_types2[un._i] != HET:
                continue

            is_phased = gt_phases1[un._i] and gt_phases2[un._i]
            # unaffected has the candidate pair on the same chromosome
            if is_phased and gt_nums1[un._i] == gt_nums2[un._i]:
                continue

            if is_phased:
                # found an unaffected with the same het-pair.
                ret['unaffected_phased'].append(un)
                if not allow_unaffected:
                    ret['candidate'] = False
            else:
                ret['unaffected_unphased'].append(un)
        if not 'candidate' in ret:
            ret['candidate'] = False
            ret['priority'] = None
        elif ret['candidate']:

            ret['priority'] = 2
            if len(ret['affected_phased']) and len(ret['unaffected_unphased']) == 0:
                ret['priority'] = 1
        return ret

    def comp_het(self, min_depth=0, gt_ll=False,
                 only_affected=True,
                 pattern_only=False):

        if pattern_only:
            af, un = empty, empty
            for kid in self.samples_with_parent:
                af |= (kid.gt_types == HET) & (kid.mom.gt_types != HOM_ALT) & (kid.dad.gt_types != HOM_ALT) \
                                            & (kid.mom.gt_types != UNKNOWN) & (kid.dad.gt_types != UNKNOWN)
        else:
            # all affecteds must be het at both sites
            af = reduce(op.or_, [s.gt_types == HET for s in self.affecteds], empty)
            # no unaffected can be homozygous alt at either site.
            un = reduce(op.and_, [s.gt_types != HOM_ALT for s in self.unaffecteds], empty)
            for kid in self.samples_with_parent:
                if not kid.affected: continue
                un &= (kid.mom.gt_types != UNKNOWN)
                un &= (kid.dad.gt_types != UNKNOWN)

            #af &= reduce(op.or_, [s.gt_types == HET for s in self.unknown], empty)

            if gt_ll:
                af &= reduce(op.and_, [s.gt_phred_ll_het <= gt_ll for s in self.affecteds])
                un &= reduce(op.and_, [s.gt_phred_ll_homalt > gt_ll for s in self.unaffecteds])

        depth = self._restrict_to_min_depth(min_depth)
        res = af & un & depth
        if res is empty:
            return 'False'
        return res

if __name__ == "__main__":

    import doctest
    HOM_REF, HET, UNKNOWN, HOM_ALT = "HOM_REF HET UNKNOWN HOM_ALT".split()
    sys.stderr.write(str(doctest.testmod(optionflags=
                     doctest.NORMALIZE_WHITESPACE
                     | doctest.ELLIPSIS
                     | doctest.REPORT_ONLY_FIRST_FAILURE, verbose=0)) + "\n")

    if len(sys.argv) > 1:
        mom = Sample('mom', False)
        dad = Sample('dad', False)
        kid = Sample('kid', True)
        kid.mom = mom
        kid.dad = dad
        fam = Family([mom, dad, kid], 'fam')

        me = fam.mendel_violations()
        for k in me:
            print(k)
            print(me[k])
            print()


        HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)

        print(fam.auto_rec(min_depth=10), "\n")

        print("auto dom:")
        print(fam.auto_dom())
        print(fam.auto_dom(min_depth=10))

        import sqlite3
        db = sqlite3.connect('test.auto_rec.db')
        fams_ped = Family.from_ped('test.auto_rec.ped')
        print(fams_ped)
