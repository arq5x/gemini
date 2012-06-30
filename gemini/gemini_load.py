#!/usr/bin/env python

# native Python imports
import os.path
import sys
import re
import collections
from optparse import OptionParser
import sqlite3
import numpy as np

# third-party imports
import cyvcf as vcf
import pysam

# gemini modules
from ped import pedformat
import infotag
import database
import annotations
import func_impact
import severe_impact
import popgen
from compression import pack_blob


def prepare_variation(args, var, v_id):
    """
    """
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
        hwe_p_value, inbreeding_coeff = popgen.get_hwe_likelihood(hom_ref, het, hom_alt, aaf)
        pi_hat = var.nucl_diversity
    else:
        aaf = extract_aaf(var)

    ########################################################
    # collect annotations from pop's custom annotation files
    ########################################################
    cyto_band  = annotations.get_cyto_info(var)
    dbsnp_info = annotations.get_dbsnp_info(var)
    in_dbsnp   = 0 if dbsnp_info.rs_ids is None else 1
    rmsk_hits  = annotations.get_rmsk_info(var)
    in_cpg     = annotations.get_cpg_island_info(var)
    in_segdup  = annotations.get_segdup_info(var)
    esp_info   = annotations.get_esp_info(var)
    # impact is a list of impacts for this variant
    impacts = None
    severe_impacts = None
    affected_gene = transcript = exon = codon_change = aa_change = consequence = effect_severity = None
    polyphen_pred = polyphen_score = sift_pred = sift_score = condel_pred = condel_score = anno_id = None
    
    if args.anno_type is not None:
        impacts = func_impact.interpret_impact(args, var)
        severe_impacts = severe_impact.interpret_severe_impact(args, var)
        if severe_impacts is not None:
            affected_gene = severe_impacts.gene
            transcript = severe_impacts.transcript
            exon = severe_impacts.exon
            codon_change = severe_impacts.codon_change
            aa_change = severe_impacts.aa_change
            consequence = severe_impacts.consequence
            effect_severity = severe_impacts.effect_severity
            polyphen_pred = severe_impacts.polyphen_pred
            polyphen_score = severe_impacts.polyphen_score
            sift_pred = severe_impacts.sift_pred
            sift_score = severe_impacts.sift_score
            condel_pred = severe_impacts.condel_pred
            condel_score = severe_impacts.condel_score
            anno_id = severe_impacts.anno_id
        #print v_id, anno_id, var.start, var.end, consequence, transcript, aa_change, affected_gene, exon, polyphen_pred
        
    # construct the filter string
    filter = None
    if var.FILTER is not None and var.FILTER != ".":
        if isinstance(var.FILTER, list):
            filter = ";".join(var.FILTER)
        else:
            filter = var.FILTER

    # build up numpy arrays for the genotype information.
    # these arrays will be pickled-to-binary, compressed,
    # and loaded as SqlLite BLOB values (see compression.pack_blob)
    gt_bases  = np.array(var.gt_bases, np.str)  # 'A/G', './.'
    gt_types  = np.array(var.gt_types, np.int8) # -1, 0, 1, 2
    gt_phases = np.array(var.gt_phases, np.bool) # T F F
    
    # were functional impacts predicted by SnpEFF or VEP?
    # if so, build up a row for each of the impacts / transcript
    variant_impacts = []
    is_exonic = False
    is_coding = False
    is_lof    = False
    gene      = None
    
    if impacts is not None:
        for idx, impact in enumerate(impacts):
            var_impact = [v_id, (idx+1), impact.gene, 
                          impact.transcript, impact.exonic, impact.coding,
                          impact.is_lof, impact.exon, impact.codon_change, 
                          impact.aa_change, impact.consequence, impact.effect_severity,
                          impact.polyphen_pred, impact.polyphen_score,
                          impact.sift_pred, impact.sift_score,
                          impact.condel_pred, impact.condel_score]
            variant_impacts.append(var_impact)
            gene = impact.gene
            if impact.exonic == True: is_exonic = True
            if impact.coding == True: is_coding = True
            if impact.is_lof == True: is_lof    = True
            
            
    # construct the core variant record.
    # 1 row per variant to VARIANTS table
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    variant = [chrom, var.start, var.end, 
               v_id, anno_id, var.REF, ','.join(var.ALT), 
               var.QUAL, filter, var.var_type, 
               var.var_subtype, pack_blob(gt_bases), pack_blob(gt_types),
               pack_blob(gt_phases), call_rate, in_dbsnp,
               dbsnp_info.rs_ids, dbsnp_info.in_omim, dbsnp_info.clin_sig,
               cyto_band, rmsk_hits, in_cpg,
               in_segdup, hom_ref, het,
               hom_alt, unknown, aaf,
               hwe_p_value, inbreeding_coeff, pi_hat,
               gene, affected_gene, transcript,    
               is_exonic, is_coding, is_lof, exon, codon_change,
               aa_change, consequence, effect_severity,
               polyphen_pred, polyphen_score, sift_pred, 
               sift_score, condel_pred, condel_score,
               infotag.get_depth(var), infotag.get_strand_bias(var), infotag.get_rms_map_qual(var),
               infotag.get_homopol_run(var), infotag.get_map_qual_zero(var), infotag.get_num_of_alleles(var),
               infotag.get_frac_dels(var), infotag.get_haplotype_score(var), infotag.get_quality_by_depth(var),
               infotag.get_allele_count(var), infotag.get_allele_bal(var), esp_info.aaf_EA, esp_info.aaf_AA,
               esp_info.aaf_ALL, esp_info.exome_chip]
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
    """
    """
    # collect of the the add'l annotation files
    annotations.load_annos()
    # open the VCF file for reading
    vcf_reader = None
    if args.vcf.endswith(".gz"):
        vcf_reader = vcf.VCFReader(open(args.vcf), 'rb', compressed=True)
    else:
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
    buffer_count = 0
    num_samples = len(vcf_reader.samples)
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

        buffer_count += 1
        # buffer full - time to insert into DB
        if buffer_count >= buffer_size:
            sys.stderr.write(str(v_id) + " variants processed.\n")
            database.insert_variation(cursor, var_buffer)
            database.insert_variation_impacts(cursor, var_impacts_buffer)
            # reset for the next batch
            var_buffer = []
            var_impacts_buffer = []
            buffer_count = 0
        v_id += 1
    # final load to the database
    database.insert_variation(cursor, var_buffer)
    database.insert_variation_impacts(cursor, var_impacts_buffer)
    sys.stderr.write(str(v_id) + " variants processed.\n")


def load(parser, args):
    if (args.db is None or args.vcf is None):
        parser.print_help()
        exit("ERROR: load needs both a VCF file and a database file\n")
    if args.anno_type not in ['snpEff', 'VEP', None]:
        parser.print_help()
        exit("\nERROR: Unsupported selection for -t\n")

    # open up a new database
    if os.path.exists(args.db):
        os.remove(args.db)
    conn = sqlite3.connect(args.db)
    conn.isolation_level = None
    c = conn.cursor()
    c.execute('PRAGMA synchronous = OFF')
    c.execute('PRAGMA journal_mode=MEMORY')
    # Create the database schema and tables.
    database.create_tables(c)
    # populate the tables.
    populate_db_from_vcf(args, c)
    # index our tables for speed
    database.create_indices(c)
    # commit data and close up
    database.close_and_commit(c, conn)
