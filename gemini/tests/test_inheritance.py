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
>>> fam.auto_rec()
True

# attach granparents to mom
>>> fam.subjects[1].mom = fam.subjects[4]
>>> fam.subjects[1].dad = fam.subjects[5]
>>> fam.auto_rec()
True
>>> _ = fam.dot()

# if grandpa is affected it is no longer autosomal recessive
>>> fam.subjects[5].affected = True
>>> fam.auto_rec()
False

>>> _ = fam.dot()

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


"""
import os
import sys
import tempfile
import atexit

from gemini import family

HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)

def tmp(pedstr, suf=".ped"):
    t = tempfile.mktemp(suffix=suf)
    atexit.register(os.unlink, t)
    with open(t, "w") as fh:
        for line in pedstr.split("\n"):
            if not line.strip(): continue
            print >> fh, line.strip()
    return t


class TestFamily(object):

    __slots__ = ('ped', 'family', 'gt_types', '_gt_types', 'gt_depths',
                 '_gt_depths', 'strict', 'subjects')

    def draw(self):
        from IPython.display import Image, display
        img = self.dot()
        return display(Image(filename=img))

    def __init__(self, ped, fam_id=None, gt_types=None, gt_depths=None):
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

    def dot(self, comment=None, path="test.gv", view=False):
        from graphviz import Digraph
        viz = Digraph(comment=comment)
        subjects = self.family.subjects
        for i, s in enumerate(subjects):
            attrs = dict(style="filled", fontcolor="white")
            attrs["color"] = {True: 'red', False: 'green', None: 'gray'}[s.affected]
            attrs["shape"] = {'male': 'square', 'female': 'circle', None: 'octagon'}[s.gender]

            if self.gt_types:
                attrs["fillcolor"] = ["white", "gray", "blue", "black"][self.gt_types[i]]
                if attrs["fillcolor"] == "black":
                    attrs["fontcolor"] = "white"
                elif attrs["fillcolor"] == "white":
                    attrs["fontcolor"] = "black"

            if s.affected is None:
                if attrs['fillcolor'] == 'gray':
                    attrs['color'] = "black"
                attrs['style'] = "filled,dashed"

            viz.node(s.name, s.name, **attrs)
        for s in subjects:
            if s.dad is not None:
                viz.edge(s.dad.name, s.name)
            if s.mom is not None:
                viz.edge(s.mom.name, s.name)
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
        def func(**kwargs):
            if 'min_depth' in kwargs:
                assert self._gt_depths is not None
            debug = kwargs.pop('debug', False)
            flt = getattr(self.family, gt)(**kwargs)
            env = {s.sample_id: i for i, s in enumerate(self.family.subjects)}
            if debug:
                print >>sys.stderr, flt
            env['gt_types'] = self.gt_types
            env['gt_depths'] = self.gt_depths
            return eval(flt, env)
        return func

if False:
    f = TestFamily("test/test.auto_rec.ped", "1")
    f.gt_types = [HET, HET, HOM_ALT]
    f.family.subjects[0].gender = "male"
    f.family.subjects[1].gender = "female"
    f.family.subjects[2].gender = "male"
    print f.family.subjects
    print f.auto_rec(strict=False)
    print f.auto_rec(strict=True)
    gm = family.Sample("grandma", False, gender="female")
    f.family.subjects[1].mom = gm
    gp = family.Sample("grandpa", None, gender="male")
    f.family.subjects[1].dad = gp
    f.gt_types.extend([HOM_REF, HET])

    f.family.subjects.extend([gm, gp])

    print f.dot("autosomal recessive")
    #f.gt_depths = [9, 9, 9]
    #print f.auto_rec(strict=True, min_depth=8)
    #print f.auto_rec(strict=True, min_depth=18)


    #f.gt_types = [HOM_ALT, HET, HOM_ALT]
    #print f.auto_rec(strict=False)
    #print f.auto_rec(strict=True)



if __name__ == "__main__":
    import sys
    import doctest
    sys.stderr.write(str(doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE, verbose=0)) + "\n")
