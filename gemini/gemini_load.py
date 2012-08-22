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

class GeminiLoader(object):
    """
    Object for creating and populating a gemini
    database and auxillary data files.
    """
    def __init__(self, args, buffer_size = 10000):
        self.args = args

        # create the gemini database
        self._create_db()
        # create a reader for the VCF file
        self.vcf_reader = self._get_vcf_reader()
        # load sample information
        self._prepare_samples()
        self.buffer_size = buffer_size
        # initialize genotype counts for each sample
        self._init_sample_gt_counts()
        
        self._get_anno_version()

    def populate_from_vcf(self):
        """
        """
        self.v_id = 1
        self.var_buffer = []
        self.var_impacts_buffer = []
        buffer_count = 0
        num_samples = len(self.samples)

        # process and load each variant in the VCF file
        for var in self.vcf_reader:
            (variant, variant_impacts) = self._prepare_variation(var)
            # add the core variant info to the variant buffer
            self.var_buffer.append(variant)
            # add each of the impact for this variant (1 per gene/transcript)
            for var_impact in variant_impacts:
                self.var_impacts_buffer.append(var_impact)
            
            # only infer genotypes if requested
            if not self.args.noload_genotypes and not self.args.no_genotypes:
                pass

            buffer_count += 1
            # buffer full - time to insert into DB
            if buffer_count >= self.buffer_size:
                sys.stderr.write(str(self.v_id) + " variants processed.\n")
                database.insert_variation(self.c, self.var_buffer)
                database.insert_variation_impacts(self.c, self.var_impacts_buffer)
                # binary.genotypes.append(var_buffer)
                # reset for the next batch
                self.var_buffer = []
                self.var_impacts_buffer = []
                buffer_count = 0
            self.v_id += 1
        # final load to the database
        database.insert_variation(self.c, self.var_buffer)
        database.insert_variation_impacts(self.c, self.var_impacts_buffer)
        sys.stderr.write(str(self.v_id) + " variants processed.\n")

    def build_indices_and_disconnect(self):
        """
        Create the db table indicies and close up
        db connection
        """
        # index our tables for speed
        database.create_indices(self.c)
        # commit data and close up
        database.close_and_commit(self.c, self.conn)

    def _get_vcf_reader(self):
        # open the VCF file for reading
        if self.args.vcf.endswith(".gz"):
            return vcf.VCFReader(open(self.args.vcf), 'rb', compressed=True)
        else:
            return vcf.VCFReader(open(self.args.vcf), 'rb')
            
    def _get_anno_version(self):
        """
        Extract the snpEff or VEP version used
        to annotate the VCF
        """
        # default to unknown version
        self.args.version = None

        if self.args.anno_type == "snpEff":
            version_string = self.vcf_reader.metadata['SnpEffVersion']
            # e.g., "SnpEff 3.0a (build 2012-07-08), by Pablo Cingolani"
            if version_string.startswith('\"SnpEff'):
                toks = version_string.split(' ')
                # e.g., make still be something like "3.0a"
                self.args.raw_version = toks[1]
                # e.g., 3.0a -> 3
                self.args.maj_version = int(self.args.raw_version.split('.')[0])
        elif self.args.anno_type == "VEP":
            pass

    def _create_db(self):
        """
        private method to open a new DB
        and create the gemini schema.
        """
        # open up a new database
        if os.path.exists(self.args.db):
            os.remove(self.args.db)
        self.conn = sqlite3.connect(self.args.db)
        self.conn.isolation_level = None
        self.c = self.conn.cursor()
        self.c.execute('PRAGMA synchronous = OFF')
        self.c.execute('PRAGMA journal_mode=MEMORY')
        # create the gemini database tables for the new DB
        database.create_tables(self.c)
        
    def _prepare_variation(self, var):
        """
        private method to collect metrics for
        a single variant (var) in a VCF file.
        """
        # these metric require that genotypes are present in the file
        call_rate = None
        hwe_p_value = None
        pi_hat = None
        inbreeding_coeff = None
        hom_ref = het = hom_alt = unknown = None
    
        # only compute certain metrics if genoypes are available
        if not self.args.no_genotypes:
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

        ############################################################
        # collect annotations from gemini's custom annotation files
        ############################################################
        cyto_band    = annotations.get_cyto_info(var)
        dbsnp_info   = annotations.get_dbsnp_info(var)
        in_dbsnp     = 0 if dbsnp_info.rs_ids is None else 1
        rmsk_hits    = annotations.get_rmsk_info(var)
        in_cpg       = annotations.get_cpg_island_info(var)
        in_segdup    = annotations.get_segdup_info(var)
        is_conserved = annotations.get_conservation_info(var)
        esp          = annotations.get_esp_info(var)
        thousandG    = annotations.get_1000G_info(var)
        recomb_rate  = annotations.get_recomb_info(var)
        gms          = annotations.get_gms(var)
        grc          = annotations.get_grc(var)
        encode_tfbs  = annotations.get_encode_tfbs(var)

        # impact is a list of impacts for this variant
        impacts = None
        severe_impacts = None
        # impact terms initialized to None for handling unannotated vcf's
        # anno_id in variants is for the transcript with the most severe impact term
        affected_gene = transcript = exon = codon_change = aa_change = aa_length = biotype = consequence = effect_severity = None
        polyphen_pred = polyphen_score = sift_pred = sift_score = anno_id = None
    
        if self.args.anno_type is not None:
            impacts = func_impact.interpret_impact(self.args, var)
            severe_impacts = severe_impact.interpret_severe_impact(self.args, var)
            if severe_impacts:
                affected_gene = severe_impacts.gene
                transcript = severe_impacts.transcript
                exon = severe_impacts.exon
                codon_change = severe_impacts.codon_change
                aa_change = severe_impacts.aa_change
                aa_length = severe_impacts.aa_length
                biotype = severe_impacts.biotype
                consequence = severe_impacts.consequence
                effect_severity = severe_impacts.effect_severity
                polyphen_pred = severe_impacts.polyphen_pred
                polyphen_score = severe_impacts.polyphen_score
                sift_pred = severe_impacts.sift_pred
                sift_score = severe_impacts.sift_score
                anno_id = severe_impacts.anno_id 
        
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
        
        # tally the genotypes
        self._update_sample_gt_counts(gt_types)
        
        # were functional impacts predicted by SnpEFF or VEP?
        # if so, build up a row for each of the impacts / transcript
        variant_impacts = []
        is_exonic = False
        is_coding = False
        is_lof    = False
        gene      = None
    
        if impacts is not None:
            for idx, impact in enumerate(impacts):
                var_impact = [self.v_id, (idx+1), impact.gene, 
                              impact.transcript, impact.exonic, impact.coding,
                              impact.is_lof, impact.exon, impact.codon_change, impact.aa_change,
                              impact.aa_length, impact.biotype, impact.consequence, impact.effect_severity,
                              impact.polyphen_pred, impact.polyphen_score,
                              impact.sift_pred, impact.sift_score]
                variant_impacts.append(var_impact)
                gene = impact.gene
                if impact.exonic == True: is_exonic = True
                if impact.coding == True: is_coding = True
                if impact.is_lof == True: is_lof    = True

            
        # construct the core variant record.
        # 1 row per variant to VARIANTS table
        chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
        variant = [chrom, var.start, var.end, 
                   self.v_id, anno_id, var.REF, ','.join(var.ALT), 
                   var.QUAL, filter, var.var_type, 
                   var.var_subtype, pack_blob(gt_bases), pack_blob(gt_types),
                   pack_blob(gt_phases), call_rate, in_dbsnp,
                   dbsnp_info.rs_ids, dbsnp_info.in_omim, dbsnp_info.clin_sig,
                   cyto_band, rmsk_hits, in_cpg,
                   in_segdup, is_conserved, hom_ref, het,
                   hom_alt, unknown, aaf,
                   hwe_p_value, inbreeding_coeff, pi_hat,
                   recomb_rate, gene, affected_gene, transcript,    
                   is_exonic, is_coding, is_lof, exon, codon_change,
                   aa_change, aa_length, biotype, consequence, effect_severity,
                   polyphen_pred, polyphen_score, sift_pred, sift_score,
                   infotag.get_depth(var), infotag.get_strand_bias(var), infotag.get_rms_map_qual(var),
                   infotag.get_homopol_run(var), infotag.get_map_qual_zero(var), infotag.get_num_of_alleles(var),
                   infotag.get_frac_dels(var), infotag.get_haplotype_score(var), infotag.get_quality_by_depth(var),
                   infotag.get_allele_count(var), infotag.get_allele_bal(var), esp.found, esp.aaf_EA, esp.aaf_AA,
                   esp.aaf_ALL, esp.exome_chip, thousandG.found, thousandG.aaf_AMR, thousandG.aaf_ASN, thousandG.aaf_AFR, 
                   thousandG.aaf_EUR, thousandG.aaf_ALL, grc, gms.illumina, gms.solid, gms.iontorrent, encode_tfbs]
        return variant, variant_impacts

    def _prepare_samples(self):
        """
        private method to load sample information
        """
        if not self.args.no_genotypes:
            self.samples = self.vcf_reader.samples
            self.sample_to_id = {}
            for idx, sample in enumerate(self.samples):
                self.sample_to_id[sample] = idx + 1

        self.ped_hash = {}
        if self.args.ped_file is not None:
            for line in open(self.args.ped_file, 'r'): 
                field = line.strip().split("\t")
                if len(field) > 1 and not field[0].startswith("#"):
                    ped = pedformat(field)
                    self.ped_hash[ped.name] = ped

        sample_list = []
        for sample in self.samples:
            i = self.sample_to_id[sample]
            if self.ped_hash.has_key(sample):
                ped = self.ped_hash[sample] 
                sample_list = [i, sample, ped.family, ped.paternal,
                               ped.maternal, ped.sex, ped.phenotype, ped.ethnicity]
            else:
                sample_list = [i, sample, None, None, None, None, None, None]
            database.insert_sample(self.c, sample_list)
    
    def _init_sample_gt_counts(self):
        """
        Initialize a 2D array of counts for tabulating
        the count of each genotype type for eaxh sample.
        
        The first dimension is one bucket for each sample.
        The second dimension (size=4) is a count for each gt type.
           Index 0 == # of hom_ref genotypes for the sample
           Index 1 == # of het genotypes for the sample
           Index 2 == # of missing genotypes for the sample
           Index 3 == # of hom_alt genotypes for the sample       
        """
        self.sample_gt_counts = np.array(np.zeros( (len(self.samples), 4) ), 
                                         dtype='uint32')

    def _update_sample_gt_counts(self, gt_types):
        """
        Update the count of each gt type for each sample
        """
        for idx, gt_type in enumerate(gt_types):
            self.sample_gt_counts[idx][gt_type] += 1

    def store_sample_gt_counts(self):
        """
        Update the count of each gt type for each sample
        """
        self.c.execute("BEGIN TRANSACTION")
        for idx, gt_counts in enumerate(self.sample_gt_counts):
            self.c.execute("""insert into sample_genotype_counts values \
                            (?,?,?,?,?)""", 
                            [idx, 
                            int(gt_counts[0]),  # hom_ref
                            int(gt_counts[1]),  # het
                            int(gt_counts[3]),  # hom_alt
                            int(gt_counts[2])]) # missing 
        self.c.execute("END")

def load(parser, args):
    if (args.db is None or args.vcf is None):
        parser.print_help()
        exit("ERROR: load needs both a VCF file and a database file\n")
    if args.anno_type not in ['snpEff', 'VEP', None]:
        parser.print_help()
        exit("\nERROR: Unsupported selection for -t\n")

    # collect of the the add'l annotation files
    annotations.load_annos()
    
    # create a new gemini loader and populate
    # the gemini db and files from the VCF
    gemini_loader = GeminiLoader(args)
    gemini_loader.populate_from_vcf()
    gemini_loader.build_indices_and_disconnect()
    gemini_loader.store_sample_gt_counts()

