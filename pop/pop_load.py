#!/usr/bin/env python

import os.path
import sys
import re
import vcf # PyVCF
from ped import pedformat
import infotag
import effect
import stats
import database
import annotations
import pysam
import sqlite3
import cPickle
import numpy as np
import zlib
import collections
from optparse import OptionParser


def get_cyto_info(var, cytoband_handle):
    """
    Returns a comma-separated list of the chromosomal
    cytobands that a variant overlaps.
    """
    # our annotation files enforce a "chr*" naming scheme.
    # enforce that the VCF file chroms abide by this.
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    cyto_band = ''
    for hit in cytoband_handle.fetch(chrom, var.start, var.end, 
                                     parser=pysam.asBed()): 
        if len(cyto_band) > 0:
            cyto_band += "," + chrom + hit.name
        else: 
            cyto_band += chrom + hit.name
    return cyto_band if len(cyto_band) > 0 else None


def get_dbsnp_info(var, dbsnp):
    """
    Returns a suite of annotations from dbSNP
    """
    DbSnpInfo = collections.namedtuple("DbSnpInfo", "rs_ids in_omim clin_sig")
    
    chrom = var.CHROM if not var.CHROM.startswith("chr") else var.CHROM[3:]
    rs_ids  = []
    clin_sigs = []
    in_omim = 0
    for hit in dbsnp.fetch(chrom, var.start, var.end, parser=pysam.asVCF()):
        rs_ids.append(hit.id)
        # load each VCF INFO key/value pair into a DICT
        info_map = {}
        for info in hit.info.split(";"):
            if info.find("=") > 0:
                (key, value) = info.split("=")
                info_map[key] = value
        # is the variant in OMIM?
        if info_map['SAO'] == 0 and info_map['OM']:
            in_omim = 1
        # what is the clinical significance of the variant?
        if info_map.get('SCS') is not None:
            clin_sigs.append(info_map['SCS'])

    rs_str = ",".join(rs_ids) if len(rs_ids) > 0 else None
    clin_sigs_str = ",".join(clin_sigs) if len(clin_sigs) > 0 else None
    return DbSnpInfo(rs_str, in_omim, clin_sigs_str)


def get_rmsk_info(var, rmsk):
    """
    Returns a comma-separated list of annotated repeats
    that overlap a variant.  Derived from the UCSC rmsk track 
    """
    # our annotation files enforce a "chr*" naming scheme.
    # enforce that the VCF file chroms abide by this.
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    rmsk_hits = []
    for hit in rmsk.fetch(chrom, var.start, var.end, parser=pysam.asBed()):
        rmsk_hits.append(hit.name)
    return ",".join(rmsk_hits) if len(rmsk_hits) > 0 else None


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
    #compute the expected number of each genotype based on p and q
    exp_hom_ref = (raf**2)*sum
    exp_het     = (2.0*(raf*aaf))*sum
    exp_hom_alt = (aaf**2)*sum
    # get the X^2 statistcs for each genotype class.
    x2_hom_ref = ((obs_hom_ref - exp_hom_ref)**2)/exp_hom_ref if exp_hom_ref > 0 else 0
    x2_hom_alt = ((obs_hom_alt - exp_hom_alt)**2)/exp_hom_alt if exp_hom_alt > 0 else 0
    x2_het     = ((obs_het - exp_het)**2)/exp_het if exp_het > 0 else 0
    x2_statistic = x2_hom_ref + x2_hom_alt + x2_het
    # return the p-value (null hyp. is that the genotypes are in HWE)
    # 1 degree of freedom b/c 3 genotypes, 2 alleles (3-2)
    
    # estimate the inbreeding coefficient (F_hat):
    # F_hat = 1 - O_hets / E_hets
    inbreeding_coeff = (1.0 - (float(obs_het)/(float(exp_het)))) if obs_het > 0 else None
    return stats.lchisqprob(x2_statistic, 1), inbreeding_coeff


def interpret_impact(var):
    """
    Interpret the report from SnpEff to determine the impact of the variant.
    For example:
    0    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T/A|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008),
    1    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T/A|OR4F5|protein_coding|CODING|ENST00000534990|exon_1_69037_69829)
    """
    try:
        effect_strings = var.INFO["EFF"].split(",")
    except KeyError:
        return None

    impact_all = [] # a list of all the transcript impacts for this variant
    for effect_string in effect_strings:
        eff_pieces = effect.eff_search.findall(effect_string)
        for piece in eff_pieces:
            impact_string = piece[0] # the predicted inpact, which is outside the ()
            impact_detail = piece[1] # all the other information, which is inside the ()
            impact_info   = effect.effect_map[impact_string]
            impact_details = effect.EffectDetails(impact_string, impact_info.priority, impact_detail)
            impact_all.append(impact_details)
    return impact_all


def prepare_variation(args, var, v_id, annos):
    
    # these metric require that genotypes are present in the file
    call_rate = None
    hwe_p_value = None
    pi_hat = None
    inbreeding_coeff = None
    hom_ref = het = hom_alt = unknown = None
    # only compute certain metrics if genoypes are available
    if not args.no_genotypes:
        hom_ref = var.num_hom_ref
        hom_alt = var.num_hom_alt
        het = var.num_het
        unknown = var.num_unknown
        call_rate = var.call_rate
        aaf = var.aaf
        hwe_p_value, inbreeding_coeff = get_hwe_likelihood(hom_ref, het, hom_alt, aaf)
        pi_hat = var.nucl_diversity
    else:
        aaf = extract_aaf(var)
    
    # collect annotations from pop's custom annotation files
    cyto_band  = get_cyto_info(var, annos['cytoband'])
    dbsnp_info = get_dbsnp_info(var, annos['dbsnp'])
    in_dbsnp   = 0 if dbsnp_info.rs_ids is None else 1
    rmsk_hits  = get_rmsk_info(var, annos['rmsk'])

    # impact is a list of impacts for this variant
    impacts = interpret_impact(var) 

    # construct the filter string
    filter = None
    if var.FILTER is not None:
        if isinstance(var.FILTER, list):
            filter = ";".join(var.FILTER)
        else:
            filter = var.FILTER

    # pack the genotype information into compressed binary values
    # for storage as BLOB fields in the database.  These binaries
    # will be decompressed and converted when SELECTed from the DB.
    gt_types  = []
    gt_phases = []
    gt_bases  = []
    for s in var.samples:
        gt_types.append(s.gt_type) if s.gt_type is not None else gt_types.append(-1)
        gt_bases.append(s.gt_bases) if s.gt_bases is not None else gt_bases.append('./.')
        gt_phases.append(s.phased) if s.phased is not None else gt_phases.append(-1)
    packed_gt_types  = sqlite3.Binary(zlib.compress(cPickle.dumps(gt_types, cPickle.HIGHEST_PROTOCOL), 9))
    packed_gt_phases = sqlite3.Binary(zlib.compress(cPickle.dumps(gt_phases, cPickle.HIGHEST_PROTOCOL), 9))
    packed_gt_bases = sqlite3.Binary(zlib.compress(cPickle.dumps(gt_bases, cPickle.HIGHEST_PROTOCOL), 9))

    variant_list = []
    # were functional impacts predicted by SnpEFF or VEP?
    if impacts is not None:
        for idx, impact in enumerate(impacts):
            var_impact = [var.CHROM, var.start, var.end, 
                          v_id, (idx+1), var.REF, 
                          ','.join(var.ALT), var.QUAL, filter,
                          var.var_type, var.var_subtype, 
                          packed_gt_bases, packed_gt_types, packed_gt_phases,
                          call_rate,
                          in_dbsnp, dbsnp_info.rs_ids, dbsnp_info.in_omim, db_snp_info.clin_sig,
                          cyto_band, rmsk_hits, 
                          hom_ref, het, 
                          hom_alt, unknown, aaf,
                          hwe_p_value, pi_hat, inbreeding_coeff,
                          impact.gene, 
                          impact.transcript, impact.exonic, impact.exon,
                          impact.coding, impact.codon_change, impact.aa_change,
                          impact.effect_name, impact.effect_severity, impact.is_lof,
                          infotag.get_depth(var), infotag.get_strand_bias(var), infotag.get_rms_map_qual(var),
                          infotag.get_homopol_run(var), infotag.get_map_qual_zero(var), infotag.get_num_of_alleles(var),
                          infotag.get_frac_dels(var), infotag.get_haplotype_score(var), infotag.get_quality_by_depth(var),
                          infotag.get_allele_count(var), infotag.get_allele_bal(var)]
            variant_list.append(var_impact)
        return variant_list
    else:
        return [[var.CHROM, var.start, var.end, 
               v_id, 1, var.REF, 
               ','.join(var.ALT), var.QUAL, filter,
               var.var_type, var.var_subtype,
               packed_gt_bases, packed_gt_types, packed_gt_phases,
               call_rate,
               in_dbsnp, dbsnp_info.rs_ids, dbsnp_info.in_omim, dbsnp_info.clin_sig,
               cyto_band, rmsk_hits, 
               hom_ref, het, 
               hom_alt, unknown, aaf,
               hwe_p_value, pi_hat, inbreeding_coeff,
               None, 
               None, None, None,
               None, None, None,
               None, None, None,
               infotag.get_depth(var), infotag.get_strand_bias(var), infotag.get_rms_map_qual(var),
               infotag.get_homopol_run(var), infotag.get_map_qual_zero(var), infotag.get_num_of_alleles(var),
               infotag.get_frac_dels(var), infotag.get_haplotype_score(var), infotag.get_quality_by_depth(var),
               infotag.get_allele_count(var), infotag.get_allele_bal(var)]]


def prepare_samples(samples, ped_file, sample_to_id, cursor):
    """
    """
    ped_hash = {}
    if ped_file is not None:
        for line in open(ped_file, 'r'): 
            field = line.strip().split("\t")
            if len(field) > 1 and not field[0].startswith("#"):
                ped = pedformat(field)
                ped_hash[ped.name] = ped

    sample_list = []
    for sample in samples:
        i = sample_to_id[sample]
        if ped_hash.has_key(sample):
            ped = ped_hash[sample] 
            sample_list = [i, sample, ped.family, ped.paternal,
                           ped.maternal, ped.sex, ped.phenotype, ped.ethnicity]
        else:
            sample_list = [i, sample, None, None, None, None, None, None]
        database.insert_sample(cursor, sample_list)


def populate_db_from_vcf(args, cursor, buffer_size = 10000):
    # collect of the the add'l annotation files
    annos = annotations.load_annos()
    # open the VCF file for reading
    vcf_reader = vcf.VCFReader(open(args.vcf), 'rb')

    if not args.no_genotypes:
        samples = vcf_reader.samples
        sample_to_id = {}
        for idx, sample in enumerate(samples):
            sample_to_id[sample] = idx + 1
        prepare_samples(samples, args.ped_file, sample_to_id, cursor)

    # load the VCF file into the variant and genotype tables
    v_id = 1
    var_buffer = []
    total_loaded = 0
    for var in vcf_reader:
        print var
        # process add'l attributes for this variant and add it to the buffer
        variant_effects = prepare_variation(args, var, v_id, annos)
        # add the impact of this variant on each gene/transcript
        for var_effect in variant_effects:
            var_buffer.append(var_effect)
        # only infer genotypes if requested
        if not args.noload_genotypes and not args.no_genotypes:
            pass
        # load the buffer to the database if it is full
        if len(var_buffer) >= buffer_size:
            total_loaded += len(var_buffer)
            sys.stderr.write(str(v_id) + " variants processed.\n")
            # add the buffers of records to the db
            database.insert_variation(cursor, var_buffer)
            # reset for the next batch
            var_buffer = []
        v_id += 1
    # final load to the database
    database.insert_variation(cursor, var_buffer)
    sys.stderr.write(str(v_id) + " variants processed.\n")


def load(parser, args):
    if (args.db is None or args.vcf is None):
        parser.print_help()
        exit()
    # open up a new database
    if os.path.exists(args.db):
        os.remove(args.db)
    conn = sqlite3.connect(args.db)
    conn.isolation_level = None
    c = conn.cursor()
    # Create the database schema and tables.
    database.create_tables(c)
    # populate the tables.
    populate_db_from_vcf(args, c)
    # index our tables for speed
    database.create_indices(c)
    # commit data and close up
    database.close_and_commit(c, conn)


if __name__ == "__main__":
    main()
