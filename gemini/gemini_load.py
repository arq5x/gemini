#!/usr/bin/env python


# native Python imports
import os.path
import sys
import sqlite3
import numpy as np

# third-party imports
import cyvcf as vcf

# gemini modules
from ped import pedformat
import infotag
import database
import annotations
import func_impact
import severe_impact
import popgen
from gemini_constants import *
from compression import pack_blob
import subprocess
from cluster_helper.cluster import cluster_view
import uuid

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

    def store_resources(self):
        """Create table of annotation resources used in this gemini database.
        """
        database.insert_resources(self.c, annotations.get_resources())

    def populate_from_vcf(self):
        """
        """
        self.v_id = 1
        self.var_buffer = []
        self.var_impacts_buffer = []
        buffer_count = 0

        # process and load each variant in the VCF file
        for var in self.vcf_reader:
            (variant, variant_impacts) = self._prepare_variation(var)
            # add the core variant info to the variant buffer
            self.var_buffer.append(variant)
            # add each of the impact for this variant (1 per gene/transcript)
            for var_impact in variant_impacts:
                self.var_impacts_buffer.append(var_impact)


            buffer_count += 1
            # buffer full - time to insert into DB
            if buffer_count >= self.buffer_size:
                sys.stderr.write(str(self.v_id) + " variants processed.\n")
                database.insert_variation(self.c, self.var_buffer)
                database.insert_variation_impacts(self.c, \
                                                  self.var_impacts_buffer)
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
        cyto_band         = annotations.get_cyto_info(var)
        rs_ids            = annotations.get_dbsnp_info(var)
        clinvar_info      = annotations.get_clinvar_info(var)
        in_dbsnp          = 0 if rs_ids is None else 1
        rmsk_hits         = annotations.get_rmsk_info(var)
        in_cpg            = annotations.get_cpg_island_info(var)
        in_segdup         = annotations.get_segdup_info(var)
        is_conserved      = annotations.get_conservation_info(var)
        esp               = annotations.get_esp_info(var)
        thousandG         = annotations.get_1000G_info(var)
        recomb_rate       = annotations.get_recomb_info(var)
        gms               = annotations.get_gms(var)
        grc               = annotations.get_grc(var)
        encode_tfbs       = annotations.get_encode_tfbs(var)
        encode_dnaseI     = annotations.get_encode_dnase_clusters(var)
        encode_cons_seg   = annotations.get_encode_consensus_segs(var)

        # impact is a list of impacts for this variant
        impacts = None
        severe_impacts = None
        # impact terms initialized to None for handling unannotated vcf's
        # anno_id in variants is for the trans. with the most severe impact term
        gene = transcript = exon = codon_change = aa_change = aa_length = \
        biotype = consequence = effect_severity = None
        is_coding = is_exonic = is_lof = 0
        polyphen_pred = polyphen_score = sift_pred = sift_score = anno_id = None

        if self.args.anno_type is not None:
            impacts = func_impact.interpret_impact(self.args, var)
            severe_impacts = \
                severe_impact.interpret_severe_impact(self.args, var)
            if severe_impacts:
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
                is_lof    = severe_impacts.is_lof

        # construct the filter string
        filter = None
        if var.FILTER is not None and var.FILTER != ".":
            if isinstance(var.FILTER, list):
                filter = ";".join(var.FILTER)
            else:
                filter = var.FILTER
        elif var.FILTER is None:
            filter = "PASS"

        # build up numpy arrays for the genotype information.
        # these arrays will be pickled-to-binary, compressed,
        # and loaded as SqlLite BLOB values (see compression.pack_blob)
        if not self.args.no_genotypes and not self.args.no_load_genotypes:
            gt_bases  = np.array(var.gt_bases, np.str)  # 'A/G', './.'
            gt_types  = np.array(var.gt_types, np.int8) # -1, 0, 1, 2
            gt_phases = np.array(var.gt_phases, np.bool) # T F F
            gt_depths = np.array(var.gt_depths, np.int32) # 10 37 0

            # tally the genotypes
            self._update_sample_gt_counts(gt_types)
        else:
            gt_bases = None
            gt_types = None
            gt_phases = None
            gt_depths = None

        # were functional impacts predicted by SnpEFF or VEP?
        # if so, build up a row for each of the impacts / transcript
        variant_impacts = []
        if impacts is not None:
            for idx, impact in enumerate(impacts):
                var_impact = [self.v_id, (idx+1), impact.gene,
                              impact.transcript, impact.is_exonic,
                              impact.is_coding, impact.is_lof,
                              impact.exon, impact.codon_change,
                              impact.aa_change, impact.aa_length,
                              impact.biotype, impact.consequence,
                              impact.effect_severity, impact.polyphen_pred,
                              impact.polyphen_score, impact.sift_pred,
                              impact.sift_score]
                variant_impacts.append(var_impact)

        # construct the core variant record.
        # 1 row per variant to VARIANTS table
        chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
        variant = [chrom, var.start, var.end,
                   self.v_id, anno_id, var.REF, ','.join(var.ALT),
                   var.QUAL, filter, var.var_type,
                   var.var_subtype, pack_blob(gt_bases), pack_blob(gt_types),
                   pack_blob(gt_phases), pack_blob(gt_depths),
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
                   cyto_band, rmsk_hits, in_cpg,
                   in_segdup, is_conserved, hom_ref, het,
                   hom_alt, unknown, aaf,
                   hwe_p_value, inbreeding_coeff, pi_hat,
                   recomb_rate, gene, transcript,
                   is_exonic, is_coding, is_lof, exon, codon_change,
                   aa_change, aa_length, biotype, consequence, effect_severity,
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
                   esp.aaf_AA, esp.aaf_ALL, esp.exome_chip, thousandG.found,
                   thousandG.aaf_AMR, thousandG.aaf_ASN, thousandG.aaf_AFR,
                   thousandG.aaf_EUR, thousandG.aaf_ALL, grc,
                   gms.illumina, gms.solid, gms.iontorrent,
                   encode_tfbs,
                   encode_dnaseI.cell_count,
                   encode_dnaseI.cell_list,
                   encode_cons_seg.gm12878,
                   encode_cons_seg.h1hesc,
                   encode_cons_seg.helas3,
                   encode_cons_seg.hepg2,
                   encode_cons_seg.huvec,
                   encode_cons_seg.k562,]
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
                               ped.maternal, ped.sex, ped.phenotype,
                               ped.ethnicity]
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
                            int(gt_counts[HOM_REF]),  # hom_ref
                            int(gt_counts[HET]),  # het
                            int(gt_counts[UNKNOWN]),  # hom_alt
                            int(gt_counts[HOM_ALT])]) # missing
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

    if use_scheduler(args):
        load_ipython(args)
    elif args.cores > 1:
        load_multicore(args)
    else:
        load_singlecore(args)


def load_singlecore(args):
    # create a new gemini loader and populate
    # the gemini db and files from the VCF
    gemini_loader = GeminiLoader(args)
    gemini_loader.store_resources()
    gemini_loader.populate_from_vcf()
    gemini_loader.build_indices_and_disconnect()

    if not args.no_genotypes and not args.no_load_genotypes:
        gemini_loader.store_sample_gt_counts()

def load_multicore(args):
    grabix_file = bgzip(args.vcf)
    chunks = load_chunks_multicore(grabix_file, args)
    merge_chunks_multicore(chunks, args.db)

def load_ipython(args):
    grabix_file = bgzip(args.vcf)
    with cluster_view(*get_ipython_args(args)) as view:
        chunks = load_chunks_ipython(grabix_file, args, view)
        merge_chunks_ipython(chunks, args.db, view)

def merge_chunks(chunks, db):
    cmd = get_merge_chunks_cmd(chunks, db)
    print "Merging chunks."
    subprocess.check_call(cmd, shell=True)
    cleanup_temp_db_files(chunks)
    return db

def get_merge_chunks_cmd(chunks, db):
    chunk_names = ""
    for chunk in chunks:
        chunk_names += " --chunkdb  " + chunk
    return "gemini merge_chunks {chunk_names} --db {db}".format(**locals())

def merge_chunks_ipython(chunks, db, view):
    if len(chunks) == 1:
        os.rename(chunks[0], db)
        return db
    else:
        sub_merges = get_chunks_to_merge(chunks)
        tmp_dbs = get_temp_dbs(len(sub_merges))
        view.map(merge_chunks, sub_merges, tmp_dbs)
        merge_chunks_ipython(tmp_dbs, db, view)

def merge_chunks_multicore(chunks, db):
    if len(chunks) == 1:
        os.rename(chunks[0], db)
        return db
    else:
        procs = []
        sub_merges = get_chunks_to_merge(chunks)
        tmp_dbs = get_temp_dbs(len(sub_merges))
        for sub_merge, tmp_db in zip(sub_merges, tmp_dbs):
            cmd = get_merge_chunks_cmd(sub_merge, tmp_db)
            procs.append(subprocess.Popen(cmd, shell=True))
        wait_until_finished(procs)
        cleanup_temp_db_files(chunks)
        merge_chunks_multicore(tmp_dbs, db)

def get_chunks_to_merge(chunks):
    sublist = list_to_sublists(chunks, 2)
    if len(sublist[-1]) > 1:
        return sublist
    else:
        sublist[-2].extend(sublist[-1])
        return sublist[:-1]

def list_to_sublists(l, n):
    """ convert list l to sublists of length n """
    return [l[i:i+n] for i in xrange(0, len(l), n)]

def get_temp_dbs(n):
    return [str(uuid.uuid4()) + ".db" for _ in xrange(n)]

def get_chunk_name(chunk):
    return "--chunkdb " + chunk

def load_chunks_multicore(grabix_file, args):
    cores = args.cores

    # specify the PED file if given one
    ped_file = ""
    if args.ped_file is not None:
        ped_file = "-p " + args.ped_file

    # specify the annotation type if given one
    anno_type = ""
    if args.anno_type is not None:
        anno_type = "-t " + args.anno_type

    submit_command = get_submit_command(args)
    vcf, _ = os.path.splitext(grabix_file)
    chunk_steps = get_chunk_steps(grabix_file, args)
    chunk_vcfs = []
    chunk_dbs = []
    procs = []
    for chunk_num, chunk in chunk_steps:
        start, stop = chunk
        print "Loading chunk " + str(chunk_num) + "." + ped_file
        gemini_load = gemini_pipe_load_cmd().format(**locals())
        procs.append(subprocess.Popen(submit_command.format(cmd=gemini_load),
                                      shell=True))

        chunk_vcf = vcf + ".chunk" + str(chunk_num)
        chunk_vcfs.append(chunk_vcf)
        chunk_dbs.append(chunk_vcf + ".db")

    wait_until_finished(procs)
    print "Done loading {0} variants in {1} chunks.".format(stop, chunk_num)
    return chunk_dbs

def load_chunks_ipython(grabix_file, args, view):
    # specify the PED file if given one
    ped_file = ""
    if args.ped_file is not None:
        ped_file = "-p " + args.ped_file

    # specify the annotation type if given one
    anno_type = ""
    if args.anno_type is not None:
        anno_type = "-t " + args.anno_type

    vcf, _ = os.path.splitext(grabix_file)
    chunk_steps = get_chunk_steps(grabix_file, args)
    total_chunks = len(chunk_steps)
    scheduler, queue, cores = get_ipython_args(args)
    load_args = {"ped_file": ped_file,
                 "anno_type": anno_type,
                 "vcf": vcf,
                 "grabix_file": grabix_file}
    chunk_dbs = view.map(load_chunk, chunk_steps, [load_args] * total_chunks)

    print "Done loading variants in {0} chunks.".format(total_chunks)
    return chunk_dbs

def load_chunk(chunk_step, kwargs):
    chunk_num, chunk = chunk_step
    start, stop = chunk
    args = combine_dicts(locals(), kwargs)
    gemini_load = gemini_pipe_load_cmd().format(**args)
    subprocess.check_call(gemini_load, shell=True)
    chunk_db = args["vcf"] + ".chunk" + str(chunk_num) + ".db"
    return chunk_db

def wait_until_finished(procs):
    [p.wait() for p in procs]

def cleanup_temp_db_files(chunk_dbs):
    for chunk_db in chunk_dbs:
        os.remove(chunk_db)

def gemini_pipe_load_cmd():
    grabix_cmd = "grabix grab {grabix_file} {start} {stop}"
    gemini_load_cmd = ("gemini load_chunk -v - {anno_type} {ped_file} "
                       "-o {start} {vcf}.chunk{chunk_num}.db")
    return " | ".join([grabix_cmd, gemini_load_cmd])

def get_chunk_steps(grabix_file, args):
    index_file = grabix_index(grabix_file)
    num_lines = get_num_lines(index_file)
    chunk_size = int(num_lines) / int(args.cores)
    print "Breaking {0} into {1} chunks.".format(grabix_file, args.cores)

    starts = []
    stops = []
    for chunk in range(0, int(args.cores)):
        start = (chunk * chunk_size) + 1
        stop  = start + chunk_size - 1
        # make sure the last chunk covers the remaining lines
        if chunk >= (args.cores - 1) and stop < num_lines:
            stop = num_lines
        starts.append(start)
        stops.append(stop)
    return list(enumerate(zip(starts, stops)))

def get_num_lines(index_file):
    with open(index_file) as index_handle:
        index_handle.next()
        num_lines = int(index_handle.next().strip())
    print "Loading %d variants." % (num_lines)
    return num_lines

def grabix_index(fname):
    if not which("grabix"):
        print_cmd_not_found_and_exit("grabix")
    index_file = fname + ".gbi"
    if file_exists(index_file):
        return index_file
    print "Indexing {0} with grabix.".format(fname)
    subprocess.check_call("grabix index {fname}".format(fname=fname), shell=True)
    return index_file

def bgzip(fname):
    if not which("bgzip"):
        print_cmd_not_found_and_exit("bgzip")
    if is_gz_file(fname):
        return fname
    bgzip_file = fname + ".gz"
    if file_exists(bgzip_file):
        return bgzip_file
    print "bgzipping {0} into {1}.".format(fname, fname + ".gz")
    subprocess.check_call("bgzip -c {fname} > {fname}.gz".format(fname=fname),
                          shell=True)
    return bgzip_file


def is_gz_file(fname):
    _, ext = os.path.splitext(fname)
    if ext == ".gz":
        return True
    else:
        return False

def get_submit_command(args):
    if args.lsf_queue:
        return get_lsf_command(args.lsf_queue)
    else:
        return "{cmd}"

def get_lsf_command(queue):
    return "bsub -K -q %s {cmd}" % (queue)

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    return os.path.exists(fname) and os.path.getsize(fname) > 0

def which(program):
    """ returns the path to an executable or None if it can't be found
     http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
     """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def combine_dicts(d1, d2):
    return dict(d1.items() + d2.items())

def get_ipython_args(args):
    if args.lsf_queue:
        return ("lsf", args.lsf_queue, args.cores)
    elif args.sge_queue:
        return ("sge", args.sge_queue, args.cores)
    elif args.torque_queue:
        return ("torque", args.torque_queue, args.cores)
    else:
        raise ValueError("ipython argument parsing failed for some reason.")

def print_cmd_not_found_and_exit(cmd):
    sys.stderr.write("Cannot find {cmd}, install it or put it in your "
                     "path.".format(cmd))
    exit(1)

def use_scheduler(args):
    if any([args.lsf_queue, args.sge_queue, args.torque_queue]):
        return True
    else:
        return False
