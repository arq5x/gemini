"""
Create a setup so we can easily define families. Input is a ped file to define
the pedigree and a vector indicating the genotype.
>>> fam = TestFamily(\"\"\"
... #family_id  sample_id   paternal_id maternal_id sex phenotype
... 1   dad   0   0   1  1
... 1   mom   0   0   2  1
... 1   kid   dad   mom   1  2
... 1   kid2   dad   mom   1  1
... 1   grandma 0   0     2  1
... 1   grandpa 0   0     1  1\"\"\")

>>> fam.gt_types = [HET, HET, HOM_ALT, HET, HET, HET]
>>> fam.gt_depths = [9] * 6
>>> _ = fam.dot()
<BLANKLINE>
auto_rec
--------
default                     True
strict=False                True
only_affected=False         True
both False                  True
<BLANKLINE>
auto_dom
--------
default                     False
strict=False                False
only_affected=False         False
both False                  False

>>> fam.auto_rec()
True

# attach granparents to mom
>>> fam.subjects[1].mom = fam.subjects[4]
>>> fam.subjects[1].dad = fam.subjects[5]
>>> fam.auto_rec()
True
>>> _ = fam.dot(tests=[])

# if grandpa is affected it is no longer autosomal recessive
>>> fam.subjects[5].affected = True
>>> fam.auto_rec()
False

>>> _ = fam.dot(tests=[])

# reset.
>>> fam.subjects[5].affected = False

# set both kids to HOM_ALT (including the
>>> fam.gt_types[3] = HOM_ALT
>>> fam.auto_rec(only_affected=True)
False
>>> fam.auto_rec(only_affected=False)
True


>>> fam.auto_rec(min_depth=10)
False

>>> fam.auto_dom()
False

# dad:un, mom:un, kid:aff, kid2:un, gma:un, gpa:un
>>> fam.gt_types = [HOM_REF, HOM_REF, HET, HET, HET, HET]
>>> fam.de_novo()
False

>>> fam.de_novo(only_affected=False)
True


>>> fam.gt_types = [HOM_ALT, HOM_REF, HET, HET, HET, HET]
>>> fam.de_novo()
False
>>> fam.gt_types = [HOM_ALT, HOM_ALT, HET, HET, HET, HET]
>>> fam.de_novo()
False

>>> fam.mendel_plausible_denovo()
True


>>> cfam = TestFamily(\"\"\"
... #family_id  sample_id   paternal_id maternal_id sex phenotype
... 1   dad   0   0   1  1
... 1   mom   0   0   2  1
... 1   akid   dad   mom   1  2
... 1   ukid   dad   mom   1  1
... 1   bkid   dad  mom   1  2\"\"\")

>>> gt_types1 = [HOM_REF, HET, HET, HOM_REF, HET]
>>> gt_bases1 = ["A/A", "A/T", "A/T", "A/A", "A/T"]

>>> gt_types2 = [HET, HOM_REF, HET, HOM_REF, HET]
>>> gt_bases2 = ["A/C", "A/A", "A/C", "A/A", "A/C"]

>>> cfam.gt_types = gt_types1

>>> cfam.comp_het()
True

>>> result = cfam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2)

# note that stuff got phased in-place:
#>>> gt_bases1, gt_bases2
#(['A/A', 'A/T', 'T|A', 'A/A', 'T|A'], ['G/C', 'G/G', 'G|C', 'G/G', 'G|C'])


>>> result['candidate']
True

>>> result['affected_phased']
[Sample(akid;affected;male), Sample(bkid;affected;male)]

>>> sorted(result.keys())
['affected_phased', 'affected_skipped', 'affected_unphased', 'candidate', 'candidates', 'priority', 'unaffected_phased', 'unaffected_unphased']


>>> assert result['affected_skipped'] == result['affected_unphased'] == result['unaffected_unphased'] == []


# remove as a candidate if even one of the affecteds doesn't share the het
# pair:
>>> gt_bases1[-1], gt_types1[-1] = "A/A", HOM_REF
>>> result = cfam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2, allow_unaffected=True)
>>> result['candidate']
False

# restore.
>>> gt_bases1[-1], gt_types1[-1] = "A/T", HET
>>> cfam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2, allow_unaffected=True)['candidate']
True


# a parent (unphaseable) has the same het pair so we know they will be phased to
# the same chrom in the affected kid.
>>> gt_types1 = [HET, HOM_REF, HET, HOM_REF, HET]
>>> gt_bases1 = ["A/T", "A/A", "A/T", "A/A", "A/T"]

>>> gt_types2 = [HET, HOM_REF, HET, HOM_REF, HET]
>>> gt_bases2 = ["G/C", "G/G", "G/C", "G/G", "G/C"]

>>> cfam.gt_types = gt_types1
>>> cfam.comp_het()
True


>>> result = cfam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2, allow_unaffected=True)
>>> gt_bases1, gt_bases2
(['A/T', 'A/A', 'A|T', 'A/A', 'A|T'], ['G/C', 'G/G', 'G|C', 'G/G', 'G|C'])

# NOTE how the variants are on the same chromosome (T, then C). so it's not a candidate.
>>> result['candidate']
False

>>> result['unaffected_unphased'], result['unaffected_phased'], result['candidate']
([Sample(dad;unaffected;male)], [], False)

# phase dad so he has same het pair (won't be a candidate):
>>> gt_bases1[0], gt_bases2[0] = "A|T", "G|C"
>>> result = cfam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2)
>>> result['unaffected_unphased'], result['unaffected_phased'], result['candidate']
([], [], False)



# unaffected kid has same het pair as affected.
>>> cfam = TestFamily(\"\"\"
... #family_id  sample_id   paternal_id maternal_id sex phenotype
... 1   dad   0   0   1  1
... 1   mom   0   0   2  1
... 1   akid   dad   mom   1  2
... 1   bkid   dad  mom   1  1\"\"\")

>>> gt_types1 = [HOM_REF, HET, HET, HET]
>>> gt_bases1 = ["A/A", "A/T", "A/T", "A/T"]

>>> gt_types2 = [HET, HOM_REF, HET, HET]
>>> gt_bases2 = ["G/C", "G/G", "G/C", "G/C"]
>>> cfam.gt_types = gt_types1

>>> cfam.comp_het()
True
>>> result = cfam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2)
>>> result['candidate']
False

# unaffected kid is ok, parent is hom_alt (this get's filtered in the first pass
# without even considering the pair).
>>> gt_bases1[-1], gt_types1[-1] = "A/A", HOM_REF
>>> gt_bases1[0], gt_types1[0] = "T/T", HOM_ALT
>>> cfam.comp_het()
False

###################################################3
# comp_het: pattern only
###################################################3

>>> cfam = TestFamily(\"\"\"
... #family_id  sample_id   paternal_id maternal_id sex phenotype
... 1   dad   0   0   1  -9
... 1   mom   0   0   2  -9
... 1   akid   dad   mom   1  -9
... 1   bkid   dad  mom   1  -9\"\"\")

>>> gt_types1 = [HOM_REF, HET, HET, HOM_REF]
>>> gt_bases1 = ["A/A", "A/T", "A/T", "A/A"]

>>> gt_types2 = [HET, HOM_REF, HET, HOM_REF]
>>> gt_bases2 = ["A/C", "A/A", "A/C", "A/A"]

>>> cfam.gt_types = gt_types1

>>> cfam.comp_het_pair(gt_types1, gt_bases1,
...                    gt_types2, gt_bases2)['candidate']
False

>>> cfam.comp_het_pair(gt_types1, gt_bases1,
...                    gt_types2, gt_bases2, pattern_only=True)
{'priority': 1, 'candidates': [Sample(akid;unknown;male)], 'candidate': True}


# get a higher priority with phased parents.
>>> gt_types1 = [HOM_REF, HET, HET, HOM_REF]
>>> gt_bases1 = ["A|A", "A|T", "A|T", "A|A"]

>>> gt_types2 = [HET, HOM_REF, HET, HOM_REF]
>>> gt_bases2 = ["A|C", "A|A", "A|C", "A|A"]
>>> cfam.comp_het_pair(gt_types1, gt_bases1,
...                    gt_types2, gt_bases2, pattern_only=True)
{'priority': 1, 'candidates': [Sample(akid;unknown;male)], 'candidate': True}

>>> fh = open('test/from_inheritance.vcf', 'w')
>>> cfam.to_vcf(fh)
>>> cfam.gt_types = gt_types2
>>> cfam.to_vcf(fh, header=False)
>>> fh.close()
>>> cfam.family.to_ped(open("test/from_inheritance.ped", "w"))

####################################################3
# auto_dom penetrance
####################################################3
>>> dfam = TestFamily(\"\"\"
... #family_id    individual_id    paternal_id    maternal_id    sex    phenotype
... 1    DS134791    DS134793    DS134792    1    2
... 1    DS134792    0    0    2    0
... 1    DS134793    0    0    1    0\"\"\")
>>> dfam.gt_types = [HET, HOM_REF, HET]
>>> dfam.auto_dom()
False

>>> dfam.auto_dom(strict=False)
True

"""

from __future__ import print_function

import os
import sys
import tempfile
import atexit

from gemini import family
import itertools as it

HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)

def tmp(pedstr, suf=".ped"):
    t = tempfile.mktemp(suffix=suf)
    atexit.register(os.unlink, t)
    with open(t, "w") as fh:
        for line in pedstr.split("\n"):
            if not line.strip(): continue
            print(line.strip(), file=fh)
    return t


class TestFamily(object):

    __slots__ = ('ped', 'family', 'gt_types', '_gt_types', 'gt_depths',
                 '_gt_depths', 'strict', 'subjects')

    def draw(self, tests=('auto_rec', 'auto_dom')):
        from IPython.display import Image, display
        if isinstance(tests, basestring):
            tests = (tests,)
        img = self.dot(tests=tests)
        return display(Image(filename=img))

    def __init__(self, ped, fam_id=None, gt_types=None, gt_depths=None):
        # can send in a family.
        if isinstance(ped, family.Family):
            self.family = ped
        else:
            if isinstance(ped, basestring) and len(ped.split("\n")) > 1:
                self.ped = tmp(ped)
            else:
                self.ped = ped
            self.family = family.Family.from_ped(self.ped) # always want 1 family

            if fam_id is None:
                assert len(self.family) == 1
                self.family = self.family.values()[0]
            else:
                self.family = self.family[fam_id]
        for s in self.family.subjects:
            if s.sample_id[0].isdigit(): s.sample_id = "s" + s.sample_id
        self.subjects = self.family.subjects

        self._gt_types = None
        self.gt_types = gt_types
        self._gt_depths = None
        self.gt_depths = gt_depths

    def dot(self, comment=None, path="test.gv", view=False, tests=('auto_rec', 'auto_dom')):
        from graphviz import Digraph
        viz = Digraph(comment=comment)
        subjects = self.family.subjects
        lookup = ["HOM_REF", "HET", "UNKOWN", "HOM_ALT"]
        for i, s in enumerate(subjects):

            attrs = dict(style="filled", fontcolor="white")
            attrs["fillcolor"] = {True: 'black', False: 'white', None: 'gray'}[s.affected]
            attrs["shape"] = {'male': 'square', 'female': 'circle', None: 'octagon'}[s.gender]

            if attrs["fillcolor"] == "black":
                attrs["fontcolor"] = "white"
            elif attrs["fillcolor"] == "white":
                attrs["fontcolor"] = "black"

            gt = lookup[self.gt_types[i]]
            label = s.name
            viz.node(s.name, label + "\n" + gt, **attrs)
        for s in subjects:
            if s.dad is not None:
                viz.edge(s.dad.name, s.name)
            if s.mom is not None:
                viz.edge(s.mom.name, s.name)
        for test in tests:
            res = {}
            res['default'] = getattr(self, test)()
            res['strict=False'] = getattr(self, test)(strict=False)
            res['only_affected=False'] = getattr(self, test)(only_affected=False)
            res['both False'] = getattr(self, test)(only_affected=False, strict=False)
            print("\n" + test)
            print("-" * len(test))
            for k in ("default", "strict=False", "only_affected=False", "both False"):
                print("%-20s\t%s" % (k, res[k]))

        viz._format = "png"
        return viz.render(path, view=view)

    @property
    def gt_types(self):
        return self._gt_types

    @gt_types.setter
    def gt_types(self, gt_types):
        if gt_types is not None:
            assert len(gt_types) == len(self.family)
            self._gt_types = gt_types

    @property
    def gt_depths(self):
        return self._gt_depths

    @gt_depths.setter
    def gt_depths(self, gt_depths):
        if gt_depths is not None:
            assert len(gt_depths) == len(self.family)
            self._gt_depths = gt_depths


    def __getattr__(self, gt):
        assert self._gt_types
        def func(*args, **kwargs):
            if 'min_depth' in kwargs:
                assert self._gt_depths is not None
            debug = kwargs.pop('debug', False)
            flt = getattr(self.family, gt)(*args, **kwargs)
            if gt == "comp_het_pair":
                return flt
            env = {s.sample_id: i for i, s in enumerate(self.family.subjects)}
            if debug:
                print(flt, file=sys.stderr)
            env['gt_types'] = self.gt_types
            env['gt_depths'] = self.gt_depths
            return eval(flt, env)
        return func

    def to_vcf(self, fh, var_dict=None, header=True, _POS=[100001]):
        if header:
            fh.write("##fileformat=VCFv4.1\n")
            fh.write("""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n""")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
            fh.write("\t".join(s.name for s in self.subjects) + "\n")
        if var_dict is None:
            var_dict = {}

        for k in ("ID", "QUAL", "INFO"):
            if k not in var_dict:
                var_dict[k] = "."
        var_dict["FILTER"] = "PASS"
        var_dict["FORMAT"] = "GT"
        if not "CHROM" in var_dict:
            var_dict["CHROM"] = "1"
        if not "POS" in var_dict:
            var_dict["POS"] = _POS[0]
            _POS[0] += 1

        if not "REF" in var_dict:
            var_dict["REF"] = "A"
        if not "ALT" in var_dict:
            var_dict["ALT"] = "G"

        # convert from number back to repr
        x = ["0/0", "0/1", "./.", "1/1"]
        formats = [x[t] for t in self.gt_types]

        if self.gt_depths:
            var_dict["FORMAT"] += ":DP"
            for i, d in enumerate(self.gt_depths):
                formats[i] += (":%d" % d)

        """
        if self.gt_phred_ll_homref:
            var_dict["FORMAT"] += ":PL"
            for i, (hom, het, alt) in enumerate(it.izip(self.gt_phred_ll_homref,
                                                self.gt_phred_ll_het,
                                                self.gt_phred_ll_homalt)):
                formats[i] += (":%s,%s,%s" % (hom, het, alt))
        """
        fh.write("{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t".format(**var_dict))
        fh.write("\t".join(formats) + "\n")

def main():
    f = TestFamily("test/test.auto_rec.ped", "1")
    f.gt_types = [HET, HET, HOM_ALT]
    f.family.subjects[0].gender = "male"
    f.family.subjects[1].gender = "female"
    f.family.subjects[2].gender = "male"
    print(f.family.subjects)
    print(f.auto_rec(strict=False))
    print(f.auto_rec(strict=True))
    gm = family.Sample("grandma", False, gender="female")
    f.family.subjects[1].mom = gm
    gp = family.Sample("grandpa", None, gender="male")
    f.family.subjects[1].dad = gp
    f.gt_types.extend([HOM_REF, HET])
    f.family.subjects.extend([gm, gp])

    print(f.dot("autosomal recessive"))
    #f.gt_depths = [9, 9, 9]
    #print f.auto_rec(strict=True, min_depth=8)
    #print f.auto_rec(strict=True, min_depth=18)


    #f.gt_types = [HOM_ALT, HET, HOM_ALT]
    #print f.auto_rec(strict=False)
    #print f.auto_rec(strict=True)
    f.to_vcf(open('a.vcf', 'w'))
    f.family.to_ped(open('a.ped', 'w'))



    import sys
    import doctest
    sys.stderr.write(str(doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE, verbose=0)) + "\n")

if __name__ == "__main__":
    main()
