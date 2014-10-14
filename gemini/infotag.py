"""
    Utility functions for extracting specific attributes
    from the VCF INFO field.
"""


def _safe_single_attr(x):
    """Handle cleaning of attributes after extraction:
      - '.'s to indicate None
      - tuples where we expect a single value
    """
    if x is None or x is ".":
        return None
    elif isinstance(x, (list, tuple)):
        return x[0]
    else:
        return x


def extract_aaf(var):
    """
    Extract the AAF directly from the INFO field.
    """
    return var.INFO.get('AF')


def get_ancestral_allele(var):
    """
    Return the reported ancestral allele if there is one
    """
    return _safe_single_attr(var.INFO.get('AA'))


def get_rms_bq(var):
    """
    Return the RMS base quality at this position
    """
    return var.INFO.get('BQ')


def get_cigar(var):
    """
    Return the cigar string describing how to align an
    alternate allele to the reference allele
    """
    return var.INFO.get('CIGAR')


def in_hm2(var):
    """
    Return whether the variant was part of HapMap2
    """
    return var.INFO.get('H2')


def in_hm3(var):
    """
    Return whether the variant was part of HapMap3
    """
    return var.INFO.get('H3')


def is_somatic(var):
    """
    Return whether the variant is somatically acquired
    """
    return var.INFO.get('SOMATIC')

def get_somatic_score(var):
    """
    Return the somatic mutation likelihood score
    """
    return var.INFO.get('SSC')


def get_depth(var):
    """
    Return the depth of aligned sequences for a given variant,
    or None if it isn't present in the VCF.
    """
    return _safe_single_attr(var.INFO.get('DP'))


def get_strand_bias(var):
    """
    Return the likelihood of strand bias,
    or None if it isn't present in the VCF.
    """
    return var.INFO.get('SB')


def get_rms_map_qual(var):
    """
    Return the RMS mapping quality,
    or None if it isn't present in the VCF.
    """
    return _safe_single_attr(var.INFO.get('MQ'))


def get_homopol_run(var):
    """
    Return the largest contiguous homopolymer run of the variant allele,
    or None if it isn't present in the VCF.
    """
    return var.INFO.get('HRun')


def get_map_qual_zero(var):
    """
    Return the total counts of mapping quality zero reads,
    or None if it isn't present in the VCF.
    """
    return _safe_single_attr(var.INFO.get('MQ0'))


def get_num_of_alleles(var):
    """
    return the total number of alleles in called genotypes,
    or None if it isn't present in the VCF.
    """
    return _safe_single_attr(var.INFO.get('AN'))


def get_frac_dels(var):
    """
    Returns the fraction of reads containing spanning deletions,
    or None if it isn't present in the VCF.
    """
    return var.INFO.get('Dels')


def get_haplotype_score(var):
    """
    Returns consistency of the site with two segregating haplotypes,
    or None if it isn't present in the VCF.
    """
    return var.INFO.get('HaplotypeScore')


def get_quality_by_depth(var):
    """
    Returns the quality by depth or the variant confidence,
    or None if it isn't present in the VCF.
    """
    return var.INFO.get('QD')


def get_allele_count(var):
    """
    Returns allele counts in genotypes,
    or None if it isn't present in the VCF.
    """
    return _safe_single_attr(var.INFO.get('AC'))


def get_allele_bal(var):
    """
    Returns allele balance for hets,
    or None if it isn't present in the VCF.
    """
    return _safe_single_attr(var.INFO.get('AB'))
