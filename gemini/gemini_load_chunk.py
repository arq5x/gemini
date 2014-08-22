#!/usr/bin/env python

# native Python imports
import os.path
import sys
import sqlite3
import numpy as np
from itertools import repeat
import json

# third-party imports
import cyvcf as vcf

# gemini modules
import version
from ped import default_ped_fields, load_ped_file
import gene_table
import infotag
import database
import annotations
import func_impact
import severe_impact
import popgen
from gemini_constants import *
from compression import pack_blob
from gemini.config import read_gemini_config


class GeminiLoader(object):
    """
    Object for creating and populating a gemini
    database and auxillary data files.
    """
    def __init__(self, args, buffer_size=10000):
        self.args = args

        # create the gemini database
        self._create_db()
        # create a reader for the VCF file
        self.vcf_reader = self._get_vcf_reader()
        # load sample information

        if not self.args.no_genotypes and not self.args.no_load_genotypes:
            # load the sample info from the VCF file.
            self._prepare_samples()
            # initialize genotype counts for each sample
            self._init_sample_gt_counts()
            self.num_samples = len(self.samples)
        else:
            self.num_samples = 0

        self.buffer_size = buffer_size
        self._get_anno_version()
        
        if not args.skip_gene_tables:
            self._get_gene_detailed()
            self._get_gene_summary()

        if self.args.anno_type == "VEP":
            self._effect_fields = self._get_vep_csq(self.vcf_reader)
        else:
            self._effect_fields = []

    def store_resources(self):
        """Create table of annotation resources used in this gemini database.
        """
        database.insert_resources(self.c, annotations.get_resources( self.args ))

    def store_version(self):
        """Create table documenting which gemini version was used for this db.
        """
        database.insert_version(self.c, version.__version__)

    def _get_vid(self):
        if hasattr(self.args, 'offset'):
            v_id = int(self.args.offset)
        else:
            v_id = 1
        return v_id

    def populate_from_vcf(self):
        """
        """
        import gemini_annotate  # avoid circular dependencies
        self.v_id = self._get_vid()
        self.counter = 0
        self.var_buffer = []
        self.var_impacts_buffer = []
        buffer_count = 0
        self.skipped = 0
        extra_file, extraheader_file = gemini_annotate.get_extra_files(self.args.db)
        extra_headers = {}
        with open(extra_file, "w") as extra_handle:
            # process and load each variant in the VCF file
            for var in self.vcf_reader:
                if self.args.passonly and (var.FILTER is not None and var.FILTER != "."):
                    self.skipped += 1
                    continue
                (variant, variant_impacts, extra_fields) = self._prepare_variation(var)
                if extra_fields:
                    extra_handle.write("%s\n" % json.dumps(extra_fields))
                    extra_headers = self._update_extra_headers(extra_headers, extra_fields)
                # add the core variant info to the variant buffer
                self.var_buffer.append(variant)
                # add each of the impact for this variant (1 per gene/transcript)
                for var_impact in variant_impacts:
                    self.var_impacts_buffer.append(var_impact)

                buffer_count += 1
                # buffer full - time to insert into DB
                if buffer_count >= self.buffer_size:
                    sys.stderr.write("pid " + str(os.getpid()) + ": " +
                                     str(self.counter) + " variants processed.\n")
                    database.insert_variation(self.c, self.var_buffer)
                    database.insert_variation_impacts(self.c,
                                                      self.var_impacts_buffer)
                    # binary.genotypes.append(var_buffer)
                    # reset for the next batch
                    self.var_buffer = []
                    self.var_impacts_buffer = []
                    buffer_count = 0
                self.v_id += 1
                self.counter += 1
        if extra_headers:
            with open(extraheader_file, "w") as out_handle:
                out_handle.write(json.dumps(extra_headers))
        else:
            os.remove(extra_file)
        # final load to the database
        self.v_id -= 1
        database.insert_variation(self.c, self.var_buffer)
        database.insert_variation_impacts(self.c, self.var_impacts_buffer)
        sys.stderr.write("pid " + str(os.getpid()) + ": " +
                         str(self.counter) + " variants processed.\n")
        if self.args.passonly:
            sys.stderr.write("pid " + str(os.getpid()) + ": " +
                             str(self.skipped) + " skipped due to having the "
                             "FILTER field set.\n")

    def _update_extra_headers(self, headers, cur_fields):
        """Update header information for extra fields.
        """
        for field, val in cur_fields.items():
            headers[field] = self._get_field_type(val, headers.get(field, "integer"))
        return headers

    def _get_field_type(self, val, cur_type):
        start_checking = False
        for name, check_fn in [("integer", int), ("float", float), ("text", str)]:
            if name == cur_type:
                start_checking = True
            if start_checking:
                try:
                    check_fn(val)
                    break
                except:
                    continue
        return name

    def build_indices_and_disconnect(self):
        """
        Create the db table indices and close up
        db connection
        """
        # index our tables for speed
        database.create_indices(self.c)
        # commit data and close up
        database.close_and_commit(self.c, self.conn)

    def _get_vcf_reader(self):
        # the VCF is a proper file
        if self.args.vcf != "-":
            if self.args.vcf.endswith(".gz"):
                return vcf.VCFReader(open(self.args.vcf), 'rb', compressed=True)
            else:
                return vcf.VCFReader(open(self.args.vcf), 'rb')
        # the VCF is being passed in via STDIN
        else:
            return vcf.VCFReader(sys.stdin, 'rb')

    def _get_anno_version(self):
        """
        Extract the snpEff or VEP version used
        to annotate the VCF
        """
        # default to unknown version
        self.args.version = None

        if self.args.anno_type == "snpEff":
            try:
                version_string = self.vcf_reader.metadata['SnpEffVersion']
            except KeyError:
                error = ("\nWARNING: VCF is not annotated with snpEff, check documentation at:\n"\
                "http://gemini.readthedocs.org/en/latest/content/functional_annotation.html#stepwise-installation-and-usage-of-snpeff\n")
                sys.exit(error)

            # e.g., "SnpEff 3.0a (build 2012-07-08), by Pablo Cingolani"
            # or "3.3c (build XXXX), by Pablo Cingolani"

            version_string = version_string.replace('"', '')  # No quotes

            toks = version_string.split()

            if "SnpEff" in toks[0]:
                self.args.raw_version = toks[1]  # SnpEff *version*, etc
            else:
                self.args.raw_version = toks[0]  # *version*, etc
            # e.g., 3.0a -> 3
            self.args.maj_version = int(self.args.raw_version.split('.')[0])

        elif self.args.anno_type == "VEP":
            pass

    def _get_vep_csq(self, reader):
        """
        Test whether the VCF header meets expectations for
        proper execution of VEP for use with Gemini.
        """
        required = ["Consequence"]
        expected = "Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE".upper()
        if 'CSQ' in reader.infos:
            parts = str(reader.infos["CSQ"].desc).split("Format: ")[-1].split("|")
            all_found = True
            for check in required:
                if check not in parts:
                    all_found = False
                    break
            if all_found:
                return parts
        # Did not find expected fields
        error = "\nERROR: Check gemini docs for the recommended VCF annotation with VEP"\
                "\nhttp://gemini.readthedocs.org/en/latest/content/functional_annotation.html#stepwise-installation-and-usage-of-vep"
        sys.exit(error)

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
        database.create_sample_table(self.c, self.args)

    def _prepare_variation(self, var):
        """private method to collect metrics for a single variant (var) in a VCF file.

        Extracts variant information, variant impacts and extra fields for annotation.
        """
        extra_fields = {}
        # these metric require that genotypes are present in the file
        call_rate = None
        hwe_p_value = None
        pi_hat = None
        inbreeding_coeff = None
        hom_ref = het = hom_alt = unknown = None

        # only compute certain metrics if genoypes are available
        if not self.args.no_genotypes and not self.args.no_load_genotypes:
            hom_ref = var.num_hom_ref
            hom_alt = var.num_hom_alt
            het = var.num_het
            unknown = var.num_unknown
            call_rate = var.call_rate
            aaf = var.aaf
            hwe_p_value, inbreeding_coeff = \
                popgen.get_hwe_likelihood(hom_ref, het, hom_alt, aaf)
            pi_hat = var.nucl_diversity
        else:
            aaf = infotag.extract_aaf(var)

        ############################################################
        # collect annotations from gemini's custom annotation files
        ############################################################
        pfam_domain = annotations.get_pfamA_domains(var)
        cyto_band = annotations.get_cyto_info(var)
        rs_ids = annotations.get_dbsnp_info(var)
        clinvar_info = annotations.get_clinvar_info(var)
        in_dbsnp = 0 if rs_ids is None else 1
        rmsk_hits = annotations.get_rmsk_info(var)
        in_cpg = annotations.get_cpg_island_info(var)
        in_segdup = annotations.get_segdup_info(var)
        is_conserved = annotations.get_conservation_info(var)
        esp = annotations.get_esp_info(var)
        thousandG = annotations.get_1000G_info(var)
        recomb_rate = annotations.get_recomb_info(var)
        gms = annotations.get_gms(var)
        grc = annotations.get_grc(var)
        in_cse = annotations.get_cse(var)
        encode_tfbs = annotations.get_encode_tfbs(var)
        encode_dnaseI = annotations.get_encode_dnase_clusters(var)
        encode_cons_seg = annotations.get_encode_consensus_segs(var)
        gerp_el = annotations.get_gerp_elements(var)
        vista_enhancers = annotations.get_vista_enhancers(var)
        cosmic_ids = annotations.get_cosmic_info(var)

        #load CADD scores by default
        if self.args.skip_cadd is False:
            (cadd_raw, cadd_scaled) = annotations.get_cadd_scores(var)
        else:
            (cadd_raw, cadd_scaled) = (None, None)

        # load the GERP score for this variant by default.
        gerp_bp = None
        if self.args.skip_gerp_bp is False:
            gerp_bp = annotations.get_gerp_bp(var)

        # impact is a list of impacts for this variant
        impacts = None
        severe_impacts = None
        # impact terms initialized to None for handling unannotated vcf's
        # anno_id in variants is for the trans. with the most severe impact term
        gene = transcript = exon = codon_change = aa_change = aa_length = \
            biotype = consequence = consequence_so = effect_severity = None
        is_coding = is_exonic = is_lof = None
        polyphen_pred = polyphen_score = sift_pred = sift_score = anno_id = None

        if self.args.anno_type is not None:
            impacts = func_impact.interpret_impact(self.args, var, self._effect_fields)
            severe_impacts = \
                severe_impact.interpret_severe_impact(self.args, var, self._effect_fields)
            if severe_impacts:
                extra_fields.update(severe_impacts.extra_fields)
                gene = severe_impacts.gene
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
                is_exonic = severe_impacts.is_exonic
                is_coding = severe_impacts.is_coding
                is_lof = severe_impacts.is_lof
                consequence_so = severe_impacts.so

        # construct the filter string
        filter = None
        if var.FILTER is not None and var.FILTER != ".":
            if isinstance(var.FILTER, list):
                filter = ";".join(var.FILTER)
            else:
                filter = var.FILTER

        vcf_id = None
        if var.ID is not None and var.ID != ".":
            vcf_id = var.ID

        # build up numpy arrays for the genotype information.
        # these arrays will be pickled-to-binary, compressed,
        # and loaded as SqlLite BLOB values (see compression.pack_blob)
        if not self.args.no_genotypes and not self.args.no_load_genotypes:
            gt_bases = np.array(var.gt_bases, np.str)  # 'A/G', './.'
            gt_types = np.array(var.gt_types, np.int8)  # -1, 0, 1, 2
            gt_phases = np.array(var.gt_phases, np.bool)  # T F F
            gt_depths = np.array(var.gt_depths, np.int32)  # 10 37 0
            gt_ref_depths = np.array(var.gt_ref_depths, np.int32)  # 2 21 0 -1
            gt_alt_depths = np.array(var.gt_alt_depths, np.int32)  # 8 16 0 -1
            gt_quals = np.array(var.gt_quals, np.float32)  # 10.78 22 99 -1

            # tally the genotypes
            self._update_sample_gt_counts(gt_types)
        else:
            gt_bases = None
            gt_types = None
            gt_phases = None
            gt_depths = None
            gt_ref_depths = None
            gt_alt_depths = None
            gt_quals = None

        if self.args.skip_info_string is False:
            info = var.INFO
        else:
            info = None

        # were functional impacts predicted by SnpEFF or VEP?
        # if so, build up a row for each of the impacts / transcript
        variant_impacts = []
        if impacts is not None:
            for idx, impact in enumerate(impacts):
                var_impact = [self.v_id, (idx + 1), impact.gene,
                              impact.transcript, impact.is_exonic,
                              impact.is_coding, impact.is_lof,
                              impact.exon, impact.codon_change,
                              impact.aa_change, impact.aa_length,
                              impact.biotype, impact.consequence,
                              impact.so, impact.effect_severity,
                              impact.polyphen_pred, impact.polyphen_score,
                              impact.sift_pred, impact.sift_score]
                variant_impacts.append(var_impact)

        # construct the core variant record.
        # 1 row per variant to VARIANTS table
        if extra_fields:
            extra_fields.update({"chrom": var.CHROM, "start": var.start, "end": var.end})
        chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
        variant = [chrom, var.start, var.end,
                   vcf_id, self.v_id, anno_id, var.REF, ','.join(var.ALT),
                   var.QUAL, filter, var.var_type,
                   var.var_subtype, pack_blob(gt_bases), pack_blob(gt_types),
                   pack_blob(gt_phases), pack_blob(gt_depths),
                   pack_blob(gt_ref_depths), pack_blob(gt_alt_depths),
                   pack_blob(gt_quals),
                   call_rate, in_dbsnp,
                   rs_ids,
                   clinvar_info.clinvar_in_omim,
                   clinvar_info.clinvar_sig,
                   clinvar_info.clinvar_disease_name,
                   clinvar_info.clinvar_dbsource,
                   clinvar_info.clinvar_dbsource_id,
                   clinvar_info.clinvar_origin,
                   clinvar_info.clinvar_dsdb,
                   clinvar_info.clinvar_dsdbid,
                   clinvar_info.clinvar_disease_acc,
                   clinvar_info.clinvar_in_locus_spec_db,
                   clinvar_info.clinvar_on_diag_assay,
                   pfam_domain, cyto_band, rmsk_hits, in_cpg,
                   in_segdup, is_conserved, gerp_bp, gerp_el,
                   hom_ref, het, hom_alt, unknown,
                   aaf, hwe_p_value, inbreeding_coeff, pi_hat,
                   recomb_rate, gene, transcript, is_exonic,
                   is_coding, is_lof, exon, codon_change, aa_change,
                   aa_length, biotype, consequence, consequence_so, effect_severity,
                   polyphen_pred, polyphen_score, sift_pred, sift_score,
                   infotag.get_ancestral_allele(var), infotag.get_rms_bq(var),
                   infotag.get_cigar(var),
                   infotag.get_depth(var), infotag.get_strand_bias(var),
                   infotag.get_rms_map_qual(var), infotag.get_homopol_run(var),
                   infotag.get_map_qual_zero(var),
                   infotag.get_num_of_alleles(var),
                   infotag.get_frac_dels(var),
                   infotag.get_haplotype_score(var),
                   infotag.get_quality_by_depth(var),
                   infotag.get_allele_count(var), infotag.get_allele_bal(var),
                   infotag.in_hm2(var), infotag.in_hm3(var),
                   infotag.is_somatic(var),
                   esp.found, esp.aaf_EA,
                   esp.aaf_AA, esp.aaf_ALL,
                   esp.exome_chip, thousandG.found,
                   thousandG.aaf_AMR, thousandG.aaf_ASN,
                   thousandG.aaf_AFR, thousandG.aaf_EUR,
                   thousandG.aaf_ALL, grc,
                   gms.illumina, gms.solid,
                   gms.iontorrent, in_cse,
                   encode_tfbs,
                   encode_dnaseI.cell_count,
                   encode_dnaseI.cell_list,
                   encode_cons_seg.gm12878,
                   encode_cons_seg.h1hesc,
                   encode_cons_seg.helas3,
                   encode_cons_seg.hepg2,
                   encode_cons_seg.huvec,
                   encode_cons_seg.k562,
                   vista_enhancers,
                   cosmic_ids,
                   pack_blob(info),
                   cadd_raw,
                   cadd_scaled]

        return variant, variant_impacts, extra_fields

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
            self.ped_hash = load_ped_file(self.args.ped_file)

        sample_list = []
        for sample in self.samples:
            i = self.sample_to_id[sample]
            if sample in self.ped_hash:
                fields = self.ped_hash[sample]
                sample_list = [i] + fields
            elif len(self.ped_hash) > 0:
                sys.exit("EXITING: sample %s found in the VCF but "
                                 "not in the PED file.\n" % (sample))
            else:
                # if there is no ped file given, just fill in the name and
                # sample_id and set the other required fields to None
                sample_list = [i, None, sample]
                sample_list += list(repeat(None, len(default_ped_fields) - 2))
            database.insert_sample(self.c, sample_list)
            
    def _get_gene_detailed(self):
        """
        define a gene detailed table
        """
        #unique identifier for each entry
        i = 0
        table_contents = detailed_list = []
        
        config = read_gemini_config( args = self.args )
        path_dirname = config["annotation_dir"]
        file_handle = os.path.join(path_dirname, 'detailed_gene_table_v75')
        
        for line in open(file_handle, 'r'):
            field = line.strip().split("\t")
            if not field[0].startswith("Chromosome"):
                i += 1
                table = gene_table.gene_detailed(field)
                detailed_list = [str(i),table.chrom,table.gene,table.is_hgnc,
                                 table.ensembl_gene_id,table.ensembl_trans_id, 
                                 table.biotype,table.trans_status,table.ccds_id, 
                                 table.hgnc_id,table.entrez,table.cds_length,table.protein_length, 
                                 table.transcript_start,table.transcript_end,
                                 table.strand,table.synonym,table.rvis,table.mam_phenotype]
                table_contents.append(detailed_list)
        database.insert_gene_detailed(self.c, table_contents)
        
    def _get_gene_summary(self):
        """
        define a gene summary table
        """
        #unique identifier for each entry
        i = 0
        contents = summary_list = []
        
        config = read_gemini_config( args = self.args )
        path_dirname = config["annotation_dir"]
        file = os.path.join(path_dirname, 'summary_gene_table_v75')
        
        for line in open(file, 'r'):
            col = line.strip().split("\t")
            if not col[0].startswith("Chromosome"):
                i += 1
                table = gene_table.gene_summary(col)
                # defaul cosmic census to False
                cosmic_census = 0
                summary_list = [str(i),table.chrom,table.gene,table.is_hgnc,
                                table.ensembl_gene_id,table.hgnc_id,
                                table.transcript_min_start,
                                table.transcript_max_end,table.strand,
                                table.synonym,table.rvis,table.mam_phenotype,
                                cosmic_census]
                contents.append(summary_list)
        database.insert_gene_summary(self.c, contents)

    def update_gene_table(self):
        """
        """
        gene_table.update_cosmic_census_genes(self.c, self.args)

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
        self.sample_gt_counts = np.array(np.zeros((len(self.samples), 4)),
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
                            int(gt_counts[HOM_REF]),  # hom_ref
                            int(gt_counts[HET]),  # het
                            int(gt_counts[HOM_ALT]),  # hom_alt
                            int(gt_counts[UNKNOWN])])  # missing
        self.c.execute("END")


def load(parser, args):
    if (args.db is None or args.vcf is None):
        parser.print_help()
        exit("ERROR: load needs both a VCF file and a database file\n")
    if args.anno_type not in ['snpEff', 'VEP', None]:
        parser.print_help()
        exit("\nERROR: Unsupported selection for -t\n")

    # collect of the the add'l annotation files
    annotations.load_annos( args )

    # create a new gemini loader and populate
    # the gemini db and files from the VCF
    gemini_loader = GeminiLoader(args)
    gemini_loader.store_resources()
    gemini_loader.store_version()

    gemini_loader.populate_from_vcf()
    gemini_loader.update_gene_table()
    # gemini_loader.build_indices_and_disconnect()

    if not args.no_genotypes and not args.no_load_genotypes:
        gemini_loader.store_sample_gt_counts()
