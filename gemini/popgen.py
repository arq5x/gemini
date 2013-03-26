import stats


def get_hwe_likelihood(obs_hom_ref, obs_het, obs_hom_alt, aaf):
    """
    Compute the likelihood of deviation from HWE using X^2,
    as well as the inbreeding coefficient.
    """
    # Bail out if aaf is undefined. This occurs
    # when there are multiple alternate alleles
    if aaf is None:
        return (None, None)

    # how many total genotypes?
    sum = (float(obs_hom_ref) + float(obs_het) + float(obs_hom_alt))
    # get the reference allele freq
    raf = 1.0 - float(aaf)
    # compute the expected number of each genotype based on p and q
    exp_hom_ref = (raf ** 2) * sum
    exp_het = (2.0 * (raf * aaf)) * sum
    exp_hom_alt = (aaf ** 2) * sum
    # get the X^2 statistcs for each genotype class.
    x2_hom_ref = ((obs_hom_ref - exp_hom_ref) ** 2) / \
        exp_hom_ref if exp_hom_ref > 0 else 0
    x2_hom_alt = ((obs_hom_alt - exp_hom_alt) ** 2) / \
        exp_hom_alt if exp_hom_alt > 0 else 0
    x2_het = ((obs_het - exp_het) ** 2) / exp_het if exp_het > 0 else 0
    x2_statistic = x2_hom_ref + x2_hom_alt + x2_het
    # return the p-value (null hyp. is that the genotypes are in HWE)
    # 1 degree of freedom b/c 3 genotypes, 2 alleles (3-2)
    # estimate the inbreeding coefficient (F_hat):
    # F_hat = 1 - O_hets / E_hets
    inbreeding_coeff = (
        1.0 - (float(obs_het) / (float(exp_het)))) if obs_het > 0 else None
    return stats.lchisqprob(x2_statistic, 1), inbreeding_coeff
