"""
Calculate the probability of a mendelian error given the genotype likelihoods
from a trio."""

import sys
from math import log10
import gzip

nan = float('nan')

class LowGenotypeException(Exception):
    pass

def rescale(li):
    s = float(sum(li))
    if s < 1e-40:
        raise LowGenotypeException
    return [v / s for v in li]

def mendelian_error(mother, father, child, pls=False):
    """
    Return the probability of a mendelian error given the log10 genotype
    likelihoods. A large value indicates a high probability of a mendelian
    error. Low values mean that the genotype-likelihoods indicate enough
    uncertainty that it could be a genotyping error.


    # everyone is het:
    >>> het = (-2.0, -0.1, -2.0)
    >>> mendelian_error(het, het, het)
    0.047...

    # parents are hom, child is het.
    >>> father = mother = [-0.6, -2.5, -2.5]
    >>> child = [-2.5, -0.6, -2.5]
    >>> mendelian_error(mother, father, child)
    0.987...

    # same as above, but more certainty in the called genotypes:
    >>> child[1] = 1e-6
    >>> mother[0] = father[0] = 1e-6
    >>> mendelian_error(mother, father, child)
    0.996...

    # everyone is confidently homozygous alt
    >>> child = father = mother = [-11.0, -11.0, -0.1]
    >>> mendelian_error(mother, father, child)
    7.55...e-11

    # everyone is less confidently homozygous refs:
    >>> child = father = mother = [-0.1, -2.0, -2.0]
    >>> mendelian_error(mother, father, child)
    0.071...

    mother and fater are homozygous alts
    >>> mother = father = [-3.0, -3.0, -0.1]

    # child is het
    >>> child = [-3., -0.1, -3.]
    >>> mendelian_error(mother, father, child)
    0.998...

    # but when the hom-alt call is close...
    >>> child = [-3., -0.1, -0.15]
    >>> mendelian_error(mother, father, child)
    0.53...

    # mother is hom_ref, dad is het, child is hom_alt
    >>> mother, father, child = (-0.1, -2, -2), (-2, -0.1, -2), (-2, -2, -0.1)
    >>> mendelian_error(mother, father, child)
    0.976...

    # mother is hom_ref, dad is hom_alt, child is hom_ref
    >>> mother, father, child = (-0.1, -2.5, -2.5), (-2.5, -2.5, -0.1), (-0.1, -2.5, -2.5)
    >>> mendelian_error(mother, father, child)
    0.993...

    # same, but child is hom_alt
    >>> mendelian_error(mother, father, (-5, -5, -0.01))
    0.994...

    # child should be het:
    >>> mendelian_error(mother, father, (-3, 0, -3))
    0.75...

    # NOTE: does oddish things if all have very low, equal values.
    >>> mendelian_error([-16.2, -16.2, -16.2], [-14.4, -15.0, -22.6], [-24.9, -21.2, -20.9])
    0.8629...

    >>> mendelian_error([-15.5, -15.8, -19.7], [-11.8, -9.9, -22.9], [-69.7, -55.9, -58.3])

    >>> mendelian_error([-3.4, -0, -2.9], [-0, -1.8, -23.0], [-6.7, 0.0, -10.7])
    0.742...

    >>> mendelian_error([34, 0, 29], [0, 18, 23], [67, 0, 107], pls=True)
    0.74...

    """
    if pls:
        mother = [m / -10.0 for m in mother]
        father = [f / -10.0 for f in father]
        child  = [c / -10.0 for c in child]
    try:
        M = rescale([10.**m for m in mother])
        F = rescale([10.**f for f in father])
        C = rescale([10.**c for c in child])
    except LowGenotypeException:
        return None

    # by ref, and alt, we mean hom_ref, hom_alt
    p_two_ref = M[0] * F[0]
    p_two_het = M[1] * F[1]
    p_two_alt = M[2] * F[2]


    # only 1 of the parents is ...
    p_one_ref = (M[0] + F[0])/2 - p_two_ref
    p_one_het = (M[1] + F[1])/2 - p_two_het
    p_one_alt = (M[2] + F[2])/2 - p_two_alt
    # divide by 2 because parents independent.

    # all options covered because, e.g. p(two_ref) == p(zero_alt)
    assert abs(sum((p_one_ref, p_one_het, p_one_alt, p_two_ref, p_two_het, p_two_alt)) - 1) < 1e-4, \
                abs(sum((p_one_ref, p_one_het, p_one_alt, p_two_ref, p_two_het, p_two_alt)) - 1)
    ##################
    # Non-violations #
    ##################
    # a. everyone is reference
    a = p_two_ref * C[0]
    # b. everyone is hom alt
    b = p_two_alt * C[2]
    # c. 1 het and 1 ref parent. child matches
    c = p_one_het * p_one_ref * (C[0] + C[1])
    # d. 1 het and 1 alt parent. child matches
    d = p_one_het * p_one_alt * (C[1] + C[2])
    # e. both parents hets. (child can be anything)
    e = p_two_het
    # f. one hom ref, one home alt. child is het
    f = p_one_ref * p_one_alt * C[1]
    #print a, b, c, d, e, f

    p_not_error = a + b + c + d + e + f
    return 1.0 - p_not_error

def xopen(f):
    return gzip.open(f) if f.endswith(".gz") else sys.stdin if "-" == f else open(f)

def main(fh, father, mother, child):

    for line in fh:
        if line.startswith("##"):
            print line,
            continue
        elif line.startswith("#CHROM"):
            print "##INFO=<ID=MEP,Number=1,Type=Float,Description=\"probability of mendelian error\">"
            print "##INFO=<ID=MER,Number=1,Type=Float,Description=\"log10 ratio of mendelian error\">"
            fields = line.rstrip().split("\t")
            samples = fields[9:]
            idxs = [9 + samples.index(s) for s in (father, mother, child)]
            print line,
            continue

        fields = line.rstrip().split("\t")
        samples = [fields[i].split(":") for i in idxs]

        fmt = fields[8].split(":")
        if "PL" in fmt:
            gli = fmt.index("PL")
            opls = [s[gli].split(",") for s in samples]
            gls = [[int(p)/-10. for p in pl] for pl in opls]
        else:
            gli = fmt.index("GL")
            ogls = [s[gli].split(",") for s in samples]
            gls = [[float(p) for g in gl] for gl in ogls]

        for i, gl in enumerate(gls):
            while sum(gls[i]) < -50:
                gls[i] = [p / 10. for p in gls[i]]
        p = mendelian_error(gls[0], gls[1], gls[2])

        if p == 1:
            mer = 100
        elif p == 0:
            mer = 0
        elif p is None:
            mer = None
        else:
            mer = log10(p / (1.0 - p))
        if p < 1 - 1e-5 or p is None:
            continue

        fields[7] += ";MEP=%.8g" % (nan if p is None else p)
        fields[7] += ";MER=%.8g" % (nan if p is None else mer)
        print "\t".join(fields)

def test():
    from random import randint
    sys.exit()

    def gen3():
        return [randint(-70, 1) / 10. for i in range(3)]

    ps = []
    for i in xrange(100000):
        a, b, c = gen3(), gen3(), gen3()
        ps.append(mendelian_error(a, b, c))
        if ps[-1] > 0.999999:
            print "mendelian error:", tuple(a), tuple(b), tuple(c)
        elif ps[-1] < 0.00001:
            print "expected       :", tuple(a), tuple(b), tuple(c)
    try:
        import pylab as pl
        pl.hist(ps, 50)
        pl.show()
    except ImportError:
        pass

def _main():
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        sys.exit(test())

    elif len(sys.argv) != 5:
        print __doc__
        print "\nUsage: %s some.vcf father_id mother_id child_id > new.vcf\n" % sys.argv[0]
        sys.exit()

    father, mother, child = sys.argv[2:]
    main(xopen(sys.argv[1]), father, mother, child)

if __name__ == "__main__":
    import doctest
    sys.stderr.write(str(doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE, verbose=0)) + "\n")
    _main()
