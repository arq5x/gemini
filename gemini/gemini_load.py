#!/usr/bin/env python

import os.path
import sys
import re
import collections
import vcf # PyVCF
from ped import pedformat
import infotag
import effect
import stats
import database
import pysam
import sqlite3
import cPickle
import numpy as np
import zlib
import annotations
from optparse import OptionParser


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


def prepare_variation(args, var, v_id):
    
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
    cyto_band  = annotations.get_cyto_info(var)
    dbsnp_info = annotations.get_dbsnp_info(var)
    in_dbsnp   = 0 if dbsnp_info.rs_ids is None else 1
    rmsk_hits  = annotations.get_rmsk_info(var)
    in_cpg     = annotations.get_cpg_island_info(var)
    in_segdup  = annotations.get_segdup_info(var)

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
    
    # were functional impacts predicted by SnpEFF or VEP?
    # if so, build up a row for each of the impacts / transcript
    variant_impacts = []
    is_exonic = False
    is_coding = False
    is_lof    = False
    if impacts is not None:
        for idx, impact in enumerate(impacts):
            var_impact = [v_id, (idx+1), impact.gene, 
                          impact.transcript, impact.exonic, impact.coding,
                          impact.is_lof, impact.exon, impact.codon_change, 
                          impact.aa_change, impact.effect_name, impact.effect_severity]
            variant_impacts.append(var_impact)
            if impact.exonic == True: is_exonic = True
            if impact.coding == True: is_coding = True
            if impact.is_lof == True: is_lof    = True

    # construct the core variant record.
    # 1 row per variant to VARIANTS table
    variant = [var.CHROM, var.start, var.end, 
               v_id, var.REF, ','.join(var.ALT), 
               var.QUAL, filter, var.var_type, 
               var.var_subtype, packed_gt_bases, packed_gt_types,
               packed_gt_phases, call_rate, in_dbsnp,
               dbsnp_info.rs_ids, dbsnp_info.in_omim, dbsnp_info.clin_sig,
               cyto_band, rmsk_hits, in_cpg,
               in_segdup, hom_ref, het,
               hom_alt, unknown, aaf,
               hwe_p_value, inbreeding_coeff, pi_hat,
               is_exonic, is_coding, is_lof,
               infotag.get_depth(var), infotag.get_strand_bias(var), infotag.get_rms_map_qual(var),
               infotag.get_homopol_run(var), infotag.get_map_qual_zero(var), infotag.get_num_of_alleles(var),
               infotag.get_frac_dels(var), infotag.get_haplotype_score(var), infotag.get_quality_by_depth(var),
               infotag.get_allele_count(var), infotag.get_allele_bal(var)]
    return variant, variant_impacts


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
    annotations.load_annos()
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
    var_impacts_buffer = []
    
    total_loaded = 0
    for var in vcf_reader:
        (variant, variant_impacts) = prepare_variation(args, var, v_id)
        # add the core variant info to the variant buffer
        var_buffer.append(variant)
        # add each of the impact for this variant (1 per gene/transcript)
        for var_impact in variant_impacts:
            var_impacts_buffer.append(var_impact)
        
        # only infer genotypes if requested
        if not args.noload_genotypes and not args.no_genotypes:
            pass
        
        # buffer full:
        # insert the buffer records into to the database
        if len(var_impacts_buffer) >= buffer_size:
            total_loaded += len(var_buffer)
            sys.stderr.write(str(v_id) + " variants processed.\n")
            database.insert_variation(cursor, var_buffer)
            database.insert_variation_impacts(cursor, var_impacts_buffer)
            # reset for the next batch
            var_buffer = []
            var_impacts_buffer = []
            
        v_id += 1
    # final load to the database
    database.insert_variation(cursor, var_buffer)
    database.insert_variation_impacts(cursor, var_impacts_buffer)
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
