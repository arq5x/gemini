import math
from collections import Counter, defaultdict
import numpy as np
from scipy.stats import binom, norm
from pandas import DataFrame
import sys

import GeminiQuery


def burden_by_gene(args):
    """
    calculates per sample the total genetic burden for each gene
    """
    query = ("SELECT gene from variants WHERE "
             "is_coding=1 and (impact_severity = 'HIGH' or "
             "polyphen_pred = 'probably_damaging')")
    _summarize_by_gene_and_sample(args, query)


def nonsynonymous_by_gene(args):
    """
    calculates per sample the total genetic burden for each gene
    """
    query = ("SELECT variant_id, gene from variants WHERE "
             "codon_change != 'None'")
    _summarize_by_gene_and_sample(args, query)

def get_calpha(args):
    """
    Calculate the C-alpha statistic for each gene based on the observed
    counts of variants in cases and controls.

    From Neale et al, PLoS Genetics, 2011.
    http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1001322
    """
    db = args.db
    if not (args.controls and args.cases):
        case, control = _get_case_and_control_samples(args)
    else:
        case = args.cases
        control = args.controls
    assert (case and control), ("Phenotypes not found in the database and "
                                "--cases and --controls are not set.")

    samples = control + case
    # p_0 = the fraction of samples that are cases (used for weighting)
    p_0 = float(len(case)) / float(len(samples))

    ns = _medium_or_high_impact_variants(args)
    variants_in_gene, variants = _calculate_counts(ns, samples)
    header = ["gene", "T", "c", "Z", "p_value"]
    print "\t".join(header)

    for gene in variants_in_gene:
        vig = variants_in_gene[gene]

        # m = the number of variants observed for this gene
        m = len(vig.keys())

        # m_n is the number of variants with n copies (i.e., samples with the variant)
        m_n = Counter([len(x) for x in vig.values()])

        # n_i is a list reflecting the total number of samples
        # having each variant
        n_i = [len(x) for x in vig.values()]

        # y_i is a list reflecting the total number of __cases__
        # having each variant
        y_i = [len(filter(lambda a: a in case, x)) for x in vig.values()]

        # "The C-alpha test statistic T contrasts the variance of each observed
        # count with the expected variance, assuming the binomial distribution."
        # In other words, given that we have n total samples and p_0 * n of them
        # are cases, we _expect_ the variant copies to be distributed among the
        # samples following a binomal distribution.  The T statistic contrasts
        # the observed count distributions with the expected:
        #
        # T = SUM{i=(1,m)} [(y_i - n_i*p_0)^2 - n_i*p_0(1 - p_0)]
        #
        T = _calculate_T(m, p_0, n_i, y_i)

        # Calculate the variance of T in order to normalize it
        c = _calculate_c(m_n, p_0)

        # The final test statistic, Z, id just the original test statistic divided
        # by its standard deviation. "We reject the null when Z is larger than expected
        # using a one-tailed standard normal distribution for reference.
        if c == 0:
            Z = np.NaN
        else:
            Z = T / math.sqrt(c)

        # sf is the survival function ... same as 1 - CDF.
        p_value = norm.sf(Z)
        # alternatie p-value 1 - scipy.stats.chi2.cdf(T**2/c, 1)
        print "\t".join([gene, str(T), str(c), str(Z), str(p_value)])


def _get_case_and_control_samples(args):
    query = ("SELECT * from samples")
    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(query)
    cases = []
    controls = []
    for row in gq:
        if int(row["phenotype"]) == 1:
            controls.append(row["name"])
        elif int(row["phenotype"]) == 2:
            cases.append(row["name"])
    return cases, controls


def _calculate_c(m_n, p_0):
    c = 0.0
    singleton_n = 0
    for n in m_n:
        if n < 2:
            singleton_n += n
            continue
        for u in xrange(n + 1):
            c += _C_term(u, n, p_0)
    if singleton_n >= 2:
        for u in xrange(singleton_n + 1):
            c += _C_term(u, singleton_n, p_0)
    return c


def _C_term(u, n, p_0):
    p_obs_u = binom(n, p_0).pmf(u)
    return ((u - n * p_0)**2 - n * p_0 * (1 - p_0))**2 * p_obs_u


def _calculate_T(m, p_0, n_i, y_i):
    T = 0.0
    singleton_n = 0
    singleton_y = 0
    for n, y in zip(n_i, y_i):
        if n < 2:
            singleton_n += n
            singleton_y += y
            continue
        T += _variant_T_term(p_0, n, y)
    if singleton_n >= 2:
        T += _variant_T_term(p_0, singleton_n, singleton_y)
    return T


def _variant_T_term(p_0, n_i, y_i):
    return (y_i - n_i * p_0)**2 - n_i * p_0 * (1 - p_0)


def _nonsynonymous_variants(db):
    query = ("SELECT variant_id, gene from variants WHERE "
             "codon_change != 'None'")
    gq = GeminiQuery.GeminiQuery(db)
    gq.run(query, show_variant_samples=True)
    return gq

def _medium_or_high_impact_variants(args):
    query = ("SELECT variant_id, gene from variants"
             " WHERE impact_severity != 'LOW'"
             " AND aaf >= %s"
             " AND aaf <= %s" % (str(args.min_aaf), str(args.max_aaf)))

    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(query, show_variant_samples=True)
    return gq

def _calculate_counts(gq, samples):
    variants = defaultdict(Counter)
    variants_in_gene = defaultdict(defaultdict)
    for row in gq:
        gene_name = row['gene']
        samples_with_variant = [x for x in row["variant_samples"] if
                                x in samples]
        if not gene_name or not samples_with_variant:
            continue
        variants_in_gene[gene_name].update({row['variant_id']:
                                            samples_with_variant})
        new_counts = Counter(samples_with_variant)
        del new_counts['']
        variants[gene_name] += new_counts
    return variants_in_gene, variants


def _summarize_by_gene_and_sample(args, query):
    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(query, show_variant_samples=True)
    burden = defaultdict(Counter)
    for row in gq:
        gene_name = row['gene']
        if not gene_name:
            continue
        new_counts = Counter(row["HET_samples"])
        # Counter can't do scalar multiplication
        new_counts = new_counts + Counter(row["HOM_ALT_samples"])
        new_counts = new_counts + Counter(row["HOM_ALT_samples"])

        del new_counts['']
        burden[gene_name] += new_counts

    df = DataFrame({})
    for gene_name, counts in burden.items():
        df = df.append(DataFrame(counts, columns=counts.keys(),
                                 index=[gene_name]))
    df = df.replace(np.NaN, 0)
    df.to_csv(sys.stdout, float_format="%d", sep="\t", index_label='gene')


def burden(parser, args):
    if args.nonsynonymous:
        nonsynonymous_by_gene(args)
    elif args.calpha:
        get_calpha(args)
    else:
        burden_by_gene(args)
