#!/usr/bin/env python
import pysam
import sqlite3
import os
import sys
import collections
import re
from unidecode import unidecode
from bx.bbi.bigwig_file import BigWigFile
from gemini.config import read_gemini_config

# dictionary of anno_type -> open Tabix file handles
annos = {}

def get_anno_files( args ):
    config = read_gemini_config( args = args )
    anno_dirname = config["annotation_dir"]
    # Default annotations -- always found
    annos = {
     'pfam_domain': os.path.join(anno_dirname, 'hg19.pfam.ucscgenes.bed.gz'),
     'cytoband': os.path.join(anno_dirname, 'hg19.cytoband.bed.gz'),
     'dbsnp': os.path.join(anno_dirname, 'dbsnp.138.vcf.gz'),
     'clinvar': os.path.join(anno_dirname, 'clinvar_20140303.vcf.gz'),
     'gwas': os.path.join(anno_dirname, 'hg19.gwas.bed.gz'),
     'rmsk': os.path.join(anno_dirname, 'hg19.rmsk.bed.gz'),
     'segdup': os.path.join(anno_dirname, 'hg19.segdup.bed.gz'),
     'conserved': os.path.join(anno_dirname, '29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed.gz'),
     'cpg_island': os.path.join(anno_dirname, 'hg19.CpG.bed.gz'),
     'dgv': os.path.join(anno_dirname, 'hg19.dgv.bed.gz'),
     'esp': os.path.join(anno_dirname,
                         'ESP6500SI.all.snps_indels.vcf.gz'),
     '1000g': os.path.join(anno_dirname,
                           'ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.2012Oct12.vcf.gz'),
     'recomb': os.path.join(anno_dirname,
                            'genetic_map_HapMapII_GRCh37.gz'),
     'gms': os.path.join(anno_dirname,
                         'GRCh37-gms-mappability.vcf.gz'),
     'grc': os.path.join(anno_dirname, 'GRC_patch_regions.bed.gz'),
     'cse': os.path.join(anno_dirname, "cse-hiseq-8_4-2013-02-20.bed.gz"),
     'encode_tfbs': os.path.join(anno_dirname,
                                 'wgEncodeRegTfbsClusteredV2.cell_count.20130213.bed.gz'),
     'encode_dnase1': os.path.join(anno_dirname,
                                   'stam.125cells.dnaseI.hg19.bed.gz'),
     'encode_consensus_segs': os.path.join(anno_dirname,
                                           'encode.6celltypes.consensus.bedg.gz'),
     'gerp_elements': os.path.join(anno_dirname, 'hg19.gerp.elements.bed.gz'),
     'vista_enhancers': os.path.join(anno_dirname, 'hg19.vista.enhancers.20131108.bed.gz'),
     'cosmic': os.path.join(anno_dirname, 'hg19.cosmic.v67.20131024.gz')
    }
    # optional annotations
    if os.path.exists(os.path.join(anno_dirname, 'hg19.gerp.bw')):
        annos['gerp_bp'] = os.path.join(anno_dirname, 'hg19.gerp.bw')
    if os.path.exists(os.path.join(anno_dirname, 'whole_genome_SNVs.tsv.compressed.gz')):
        annos['cadd_score'] = os.path.join(anno_dirname, 'whole_genome_SNVs.tsv.compressed.gz')
    return annos

class ClinVarInfo(object):
    def __init__(self):
        self.clinvar_dbsource = None
        self.clinvar_dbsource_id = None
        self.clinvar_origin = None
        self.clinvar_sig = None
        self.clinvar_dsdb = None
        self.clinvar_dsdbid = None
        self.clinvar_disease_name = None
        self.clinvar_disease_acc = None
        self.clinvar_in_omim = None
        self.clinvar_in_locus_spec_db = None
        self.clinvar_on_diag_assay = None

        self.origin_code_map = {'0': 'unknown',
                                '1': 'germline',
                                '2': 'somatic',
                                '4': 'inherited',
                                '8': 'paternal',
                                '16': 'maternal',
                                '32': 'de-novo',
                                '64': 'biparental',
                                '128': 'uniparental',
                                '256': 'not-tested',
                                '512': 'tested-inconclusive',
                                '1073741824': 'other'}

        self.sig_code_map = {'0': 'unknown',
                             '1': 'untested',
                             '2': 'non-pathogenic',
                             '3': 'probable-non-pathogenic',
                             '4': 'probable-pathogenic',
                             '5': 'pathogenic',
                             '6': 'drug-response',
                             '7': 'histocompatibility',
                             '255': 'other'}

    def __repr__(self):
        return '\t'.join([self.clinvar_dbsource,
                          self.clinvar_dbsource_id,
                          self.clinvar_origin,
                          self.clinvar_sig,
                          self.clinvar_dsdb,
                          self.clinvar_dsdbid,
                          self.clinvar_disease_name,
                          self.clinvar_disease_acc,
                          str(self.clinvar_in_omim),
                          str(self.clinvar_in_locus_spec_db),
                          str(self.clinvar_on_diag_assay)])

    def lookup_clinvar_origin(self, origin_code):
        try:
            return self.origin_code_map[origin_code]
        except KeyError:
            return None

    def lookup_clinvar_significance(self, sig_code):
        if "|" not in sig_code:
            try:
                return self.sig_code_map[sig_code]
            except KeyError:
                return None
        else:
            sigs = set(sig_code.split('|'))
            # e.g., 255|255|255
            if len(sigs) == 1:
                try:
                    return self.sig_code_map[sigs.pop()]
                except KeyError:
                    return None
            # e.g., 1|5|255
            else:
                return "mixed"


ESPInfo = collections.namedtuple("ESPInfo",
                                 "found \
                                  aaf_EA \
                                  aaf_AA \
                                  aaf_ALL \
                                  exome_chip")
ENCODEDnaseIClusters = collections.namedtuple("ENCODEDnaseIClusters",
                                              "cell_count \
                                         cell_list")
ENCODESegInfo = collections.namedtuple("ENCODESegInfo",
                                       "gm12878 \
                                         h1hesc \
                                         helas3 \
                                         hepg2 \
                                         huvec \
                                         k562")
ThousandGInfo = collections.namedtuple("ThousandGInfo",
                                       "found \
                                        aaf_ALL \
                                        aaf_AMR \
                                        aaf_ASN \
                                        aaf_AFR \
                                        aaf_EUR")

def load_annos( args ):
    """
    Populate a dictionary of Tabixfile handles for
    each annotation file.  Other modules can then
    access a given handle and fetch data from it
    as follows:

    dbsnp_handle = annotations.annos['dbsnp']
    hits = dbsnp_handle.fetch(chrom, start, end)
    """
    anno_files = get_anno_files( args )
    for anno in anno_files: 
        try:
            # .gz denotes Tabix files.
            if anno_files[anno].endswith(".gz"):
                annos[anno] = pysam.Tabixfile(anno_files[anno])
            # .bw denotes BigWig files.
            elif anno_files[anno].endswith(".bw"):
                annos[anno] = BigWigFile( open( anno_files[anno] ) )

        except IOError:
            sys.exit("Gemini cannot open this annotation file: %s. \n"
                     "Have you installed the annotation files?  If so, "
                     "have they been moved or deleted? Exiting...\n\n"
                     "For more details:\n\t"
                     "http://gemini.readthedocs.org/en/latest/content/"
                     "#installation.html\#installing-annotation-files\n"
                     % anno_files[anno])

# ## Standard access to Tabix indexed files


def _get_hits(coords, annotation, parser_type):
    """Retrieve BED information, recovering if BED annotation file does have a chromosome.
    """
    if parser_type == "bed":
        parser = pysam.asBed()
    elif parser_type == "vcf":
        parser = pysam.asVCF()
    elif parser_type == "tuple":
        parser = pysam.asTuple()
    elif parser_type is None:
        parser = None
    else:
        raise ValueError("Unexpected parser type: %s" % parser)
    chrom, start, end = coords
    try:
        hit_iter = annotation.fetch(str(chrom), start, end, parser=parser)
    # catch invalid region errors raised by ctabix
    except ValueError:
        hit_iter = []
    # recent versions of pysam return KeyError
    except KeyError:
        hit_iter = []
    return hit_iter

def _get_bw_summary(coords, annotation):
    """Return summary of BigWig scores in an interval
    """
    chrom, start, end = coords
    try:
        return annotation.summarize(str(chrom), start, end, end-start).min_val[0]
    except AttributeError:
        return None


def _get_chr_as_grch37(chrom):
    if chrom in ["chrM"]:
        return "MT"
    return chrom if not chrom.startswith("chr") else chrom[3:]


def _get_chr_as_ucsc(chrom):
    return chrom if chrom.startswith("chr") else "chr" + chrom


def guess_contig_naming(anno):
    """Guess which contig naming scheme a given annotation file uses.
    """
    chr_names = [x for x in anno.contigs if x.startswith("chr")]
    if len(chr_names) > 0:
        return "ucsc"
    else:
        return "grch37"


def _get_var_coords(var, naming):
    """Retrieve variant coordinates from multiple input objects.
    """
    if isinstance(var, dict) or isinstance(var, sqlite3.Row):
        chrom = var["chrom"]
        start = int(var["start"])
        end = int(var["end"])
    else:
        chrom = var.CHROM
        start = var.start
        end = var.end
    if naming == "ucsc":
        chrom = _get_chr_as_ucsc(chrom)
    elif naming == "grch37":
        chrom = _get_chr_as_grch37(chrom)
    return chrom, start, end

def _get_cadd_scores(var, labels, hit):
    """
    get cadd scores
    """
    raw = hit[3].split(",")
    scaled = hit[4].split(",")
    
    p = re.compile(str(var.ALT[0]))
    for m in p.finditer(str(labels[hit[2]])):
        pos = m.start()
        return raw[pos], scaled[pos]
        

def annotations_in_region(var, anno, parser_type=None, naming="ucsc"):
    """Iterator of annotations found in a genomic region.

    - var: PyVCF object or database query with chromosome, start and end.
    - anno: pysam Tabix annotation file or string to reference
            a standard annotation
    - parser_type: string specifying the filetype of the tabix file
    - naming: chromosome naming scheme used, ucsc or grch37
    """
    coords = _get_var_coords(var, naming)
    if isinstance(anno, basestring):
        anno = annos[anno]
    return _get_hits(coords, anno, parser_type)


def bigwig_summary(var, anno, naming="ucsc"):
    coords = _get_var_coords(var, naming)
    if isinstance(anno, basestring):
        anno = annos[anno]
    return _get_bw_summary(coords, anno)


# ## Track-specific annotations
def get_cpg_island_info(var):
    """
    Returns a boolean indicating whether or not the
    variant overlaps a CpG island
    """
    for hit in annotations_in_region(var, "cpg_island", "bed"):
        return True
    return False

# def get_dbNSFP_info(var, impacts):
#     """
#     Returns Polyphen, SIFT, etc. from dbNSFP annotation file.
#     One prediction per transcript.
    
#     LIMITATION: only handles bi-allelic loci
#     """

#     # is this variant predicted to be nonsynonymous for any of the transcripts?
#     # if not, we can skip dnNSFP.
#     non_syn_impacts = [imp for imp in impacts \
#                                if imp.consequence == 'non_syn_coding']

#     if len(non_syn_impacts) > 0:
#         for hit in annotations_in_region(var, "dbnsfp", parser_type="tuple", naming="grch37"):

#             if var.POS == int(hit[1]) and \
#                var.REF == hit[2] and \
#                var.ALT[0] == hit[3]:
                
#                 transcripts = hit[7].split(';')
#                 aapos = hit[8].split(';')
#                 pp_scores = hit[11].split(';')

#                 if len(transcripts) != len(pp_scores):
#                     print var.POS, var.REF, var.ALT[0], [i.transcript for i in non_syn_impacts], \
#                           [i.polyphen_pred for i in non_syn_impacts], [i.polyphen_score for i in non_syn_impacts], \
#                           hit[7], hit[8], hit[11], hit[12]
#     else:
#         pass

def get_cyto_info(var):
    """
    Returns a comma-separated list of the chromosomal
    cytobands that a variant overlaps.
    """
    cyto_band = ''
    for hit in annotations_in_region(var, "cytoband", "bed"):
        if len(cyto_band) > 0:
            cyto_band += "," + hit.contig + hit.name
        else:
            cyto_band += hit.contig + hit.name
    return cyto_band if len(cyto_band) > 0 else None

def get_gerp_bp(var):
    """
    Returns a summary of the GERP scores for the variant.
    """
    if "gerp_bp" not in annos:
        raise IOError("Need to download BigWig file with GERP scores per base pair. "
                      "Run `gemini update --dataonly --extra gerp_bp")
    gerp = bigwig_summary(var, "gerp_bp")
    return gerp

def get_gerp_elements(var):
    """
    Returns the GERP element information.
    """
    p_vals = []
    for hit in annotations_in_region(var, "gerp_elements", "tuple"):
        p_vals.append(hit[3])
    if len(p_vals) == 1:
        return p_vals[0]
    elif len(p_vals) > 1:
        return min(float(p) for p in p_vals)
    else:
        return None

def get_vista_enhancers(var):
    """
    Returns the VISTA enhancer information.
    """
    vista_enhancers = []
    for hit in annotations_in_region(var, "vista_enhancers", "tuple"):
        vista_enhancers.append(hit[4])
    return ",".join(vista_enhancers) if len(vista_enhancers) > 0 else None


def get_cadd_scores(var):
    """
    Returns the C-raw scores & scaled scores (CADD) to predict deleterious
    variants. Implemented only for SNV's
    """
    if "cadd_score" not in annos:
        raise IOError("Need to download the CADD data file for deleteriousness."
                      "Run `gemini update --dataonly --extra cadd_score")

    cadd_raw = cadd_scaled = None
    labels = {"A":"CGT", "C":"AGT", "G":"ACT", "T":"ACG", "R":"ACGT", "M":"ACGT"}
    
    for hit in annotations_in_region(var, "cadd_score", "tuple", "grch37"):
        # we want exact position mapping here and not a range (end-start) as
        # returned in hit (e.g. indels) & we do not want to consider del & ins
        if str(hit[1]) == str(var.POS) and len(var.REF) == 1 and \
           len(var.ALT[0]) == 1:
           
            if str(hit[2]) == var.REF and str(var.ALT[0]) in labels[hit[2]]:
               (cadd_raw, cadd_scaled) = _get_cadd_scores(var, labels, hit)
            
            # consider ref cases with ambiguity codes R (G,A) and M (A,C)
            elif ((str(hit[2]) == 'R'  and var.REF in('G','A')) or \
                (str(hit[2]) == 'M'  and var.REF in('A','C'))) and \
                str(var.ALT[0]) in labels[hit[2]]:
                (cadd_raw, cadd_scaled) = _get_cadd_scores(var, labels, hit)
    
    return (cadd_raw, cadd_scaled)
     
    
def get_pfamA_domains(var):
    """
    Returns pfamA domains that a variant overlaps
    """
    pfam_domain = []
    for hit in annotations_in_region(var, "pfam_domain", "bed"):
        pfam_domain.append(hit.name)
    return ",".join(pfam_domain) if len(pfam_domain) > 0 else None


def get_cosmic_info(var):
    """
    Returns a list of COSMIC ids associated with given variant

    E.g. from COSMIC VCF
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
    chrM    1747    COSN408408      G       A       .       .       .
    chrM    2700    COSN408409      G       A       .       .       .
    chr1    42880262    COSM464635  G   C   .   .   AA=p.D224H;CDS=c.670G>C;CNT=1;GENE=RIMKLA;STRAND=+
    chr1    42880269    COSM909628  G   A   .   .   AA=p.G226D;CDS=c.677G>A;CNT=1;GENE=RIMKLA;STRAND=+
    chr1    42880284    COSM1502979 G   T   .   .   AA=p.C231F;CDS=c.692G>T;CNT=1;GENE=RIMKLA;STRAND=+
    chr1    42880318    COSM681351  T   A   .   .   AA=p.F242L;CDS=c.726T>A;CNT=1;GENE=RIMKLA;STRAND=+
    chr1    42880337    COSM464636  G   A   .   .   AA=p.D249N;CDS=c.745G>A;CNT=1;GENE=RIMKLA;STRAND=+
    chr1    42880384    COSM909629  T   C   .   .   AA=p.N264N;CDS=c.792T>C;CNT=1;GENE=RIMKLA;STRAND=+
    chr1    42880415    COSM909630  G   C   .   .   AA=p.G275R;CDS=c.823G>C;CNT=1;GENE=RIMKLA;STRAND=+
    """
    # report the first overlapping ClinVar variant Most often, just one).
    cosmic_ids = []
    for hit in annotations_in_region(var, "cosmic", "vcf", "ucsc"):
        cosmic_ids.append(hit.id)
    return ",".join(cosmic_ids) if len(cosmic_ids) > 0 else None


def get_clinvar_info(var):
    """
    Returns a suite of annotations from ClinVar

    ClinVarInfo named_tuple:
    --------------------------------------------------------------------------
    # clinvar_dbsource         = CLNSRC=OMIM Allelic Variant;
    # clinvar_dbsource_id      = CLNSRCID=103320.0001;
    # clinvar_origin           = CLNORIGIN=1
    # clinvar_sig              = CLNSIG=5
    # clinvar_dsdb             = CLNDSDB=GeneReviews:NCBI:OMIM:Orphanet;
    # clinvar_dsdbid           = CLNDSDBID=NBK1168:C1850792:254300:590;
    # clinvar_disease_name     = CLNDBN=Myasthenia\x2c limb-girdle\x2c familial;
    # clinvar_disease_acc      = CLNACC=RCV000019902.1
    # clinvar_in_omim          = OM
    # clinvar_in_locus_spec_db = LSD
    # clinvar_on_diag_assay    = CDA
    """

    clinvar = ClinVarInfo()

    # report the first overlapping ClinVar variant Most often, just one).
    for hit in annotations_in_region(var, "clinvar", "vcf", "grch37"):
        # load each VCF INFO key/value pair into a DICT
        info_map = {}
        for info in hit.info.split(";"):
            if info.find("=") > 0:
                (key, value) = info.split("=")
                info_map[key] = value
            else:
                info_map[info] = True

        raw_dbsource = info_map['CLNSRC'] or None
        #interpret 8-bit strings and convert to plain text
        clinvar.clinvar_dbsource = unidecode(raw_dbsource.decode('utf-8'))
        clinvar.clinvar_dbsource_id = info_map['CLNSRCID'] or None
        clinvar.clinvar_origin           = \
            clinvar.lookup_clinvar_origin(info_map['CLNORIGIN'])
        clinvar.clinvar_sig              = \
            clinvar.lookup_clinvar_significance(info_map['CLNSIG'])
        clinvar.clinvar_dsdb = info_map['CLNDSDB'] or None
        clinvar.clinvar_dsdbid = info_map['CLNDSDBID'] or None
        # Remap all unicode characters into plain text string replacements
        raw_disease_name = info_map['CLNDBN'] or None
        clinvar.clinvar_disease_name = unidecode(raw_disease_name.decode('utf-8'))
        # Clinvar represents commas as \x2c.  Make them commas.
        clinvar.clinvar_disease_name = clinvar.clinvar_disease_name.decode('string_escape')

        clinvar.clinvar_disease_acc = info_map['CLNACC'] or None
        clinvar.clinvar_in_omim = 1 if 'OM' in info_map else 0
        clinvar.clinvar_in_locus_spec_db = 1 if 'LSD' in info_map else 0
        clinvar.clinvar_on_diag_assay = 1 if 'CDA' in info_map else 0

    return clinvar


def get_dbsnp_info(var):
    """
    Returns a suite of annotations from dbSNP
    """
    rs_ids = []
    for hit in annotations_in_region(var, "dbsnp", "vcf", "grch37"):
        rs_ids.append(hit.id)
        # load each VCF INFO key/value pair into a DICT
        info_map = {}
        for info in hit.info.split(";"):
            if info.find("=") > 0:
                (key, value) = info.split("=")
                info_map[key] = value

    return ",".join(rs_ids) if len(rs_ids) > 0 else None


def get_esp_info(var):
    """
    Returns a suite of annotations from the ESP project

    ESP reports the minor allele frequency (MAF), not the
    alternate allele frequency (AAF). We must therefore figure
    out whther the reference or alternate allele is the minor allele.

    1       69496   rs150690004     G       A       .       PASS    DBSNP=dbSNP_134;EA_AC=2,6764;AA_AC=23,3785;TAC=25,10549;MAF=0.0296,0.604,0.2364;GTS=AA,AG,GG;EA_GTC=0,2,3381;AA_GTC=5,13,1886;GTC=5,15,5267;DP=91;GL=OR4F5;CP=0.5;CG=2.3;AA=G;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;GM=NM_001005484.1;FG=missense;AAC=SER/GLY;PP=136/306;CDP=406;GS=56;PH=benign
    1       69511   rs75062661      A       G       .       PASS    DBSNP=dbSNP_131;EA_AC=5337,677;AA_AC=1937,1623;TAC=7274,2300;MAF=11.2571,45.5899,24.0234;GTS=GG,GA,AA;EA_GTC=2430,477,100;AA_GTC=784,369,627;GTC=3214,846,727;DP=69;GL=OR4F5;CP=1.0;CG=1.1;AA=G;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;GM=NM_001005484.1;FG=missense;AAC=ALA/THR;PP=141/306;CDP=421;GS=58;PH=benign
    """
    aaf_EA = aaf_AA = aaf_ALL = None
    maf = fetched = con = []
    exome_chip = False
    found = False
    info_map = {}
    for hit in annotations_in_region(var, "esp", "vcf", "grch37"):
        if hit.contig not in ['Y']:
            fetched.append(hit)
            # We need a single ESP entry for a variant
            if fetched != None and len(fetched) == 1 and \
                    hit.alt == var.ALT[0] and hit.ref == var.REF:
                found = True
                # loads each VCF INFO key/value pair into a DICT
                for info in hit.info.split(";"):
                    if info.find("=") > 0:
                    # splits on first occurence of '='
                    # useful to handle valuerror: too many values to unpack (e.g (a,b) = split(",", (a,b,c,d)) for cases like
                    # SA=http://www.ncbi.nlm.nih.gov/sites/varvu?gene=4524&amp%3Brs=1801131|http://omim.org/entry/607093#0004
                        (key, value) = info.split("=", 1)
                        info_map[key] = value

                # get the allele counts so that we can compute alternate allele frequencies
                # example: EA_AC=2,6764;AA_AC=23,3785;TAC=25,10549
                if info_map.get('EA_AC') is not None:
                    lines = info_map['EA_AC'].split(",")
                    aaf_EA = float(lines[0]) / (float(lines[0]) + float(lines[1]))

                if info_map.get('AA_AC') is not None:
                    lines = info_map['AA_AC'].split(",")
                    aaf_AA = float(lines[0]) / (float(lines[0]) + float(lines[1]))

                if info_map.get('TAC') is not None:
                    lines = info_map['TAC'].split(",")
                    aaf_ALL = float(lines[0]) / (float(lines[0]) + float(lines[1]))

                # Is the SNP on an human exome chip?
                if info_map.get('EXOME_CHIP') is not None and \
                        info_map['EXOME_CHIP'] == "no":
                    exome_chip = 0
                elif info_map.get('EXOME_CHIP') is not None and \
                        info_map['EXOME_CHIP'] == "yes":
                    exome_chip = 1

    return ESPInfo(found, aaf_EA, aaf_AA, aaf_ALL, exome_chip)


def get_1000G_info(var):
    """
    Returns a suite of annotations from the 1000 Genomes project
    """
    #fetched = []
    info_map = {}
    found = False

    for hit in annotations_in_region(var, "1000g", "vcf", "grch37"):
        # We need to ensure we are dealing with the exact sample variant
        # based on position and the alleles present.
        if var.start == hit.pos and \
           var.ALT[0] == hit.alt and \
           hit.ref == var.REF:
            for info in hit.info.split(";"):
                if info.find("=") > 0:
                    (key, value) = info.split("=", 1)
                    info_map[key] = value
            found = True

    return ThousandGInfo(found, info_map.get('AF'), info_map.get('AMR_AF'),
                         info_map.get('ASN_AF'), info_map.get('AFR_AF'),
                         info_map.get('EUR_AF'))


def get_rmsk_info(var):
    """
    Returns a comma-separated list of annotated repeats
    that overlap a variant.  Derived from the UCSC rmsk track
    """
    rmsk_hits = []
    for hit in annotations_in_region(var, "rmsk", "bed"):
        rmsk_hits.append(hit.name)
    return ",".join(rmsk_hits) if len(rmsk_hits) > 0 else None


def get_segdup_info(var):
    """
    Returns a boolean indicating whether or not the
    variant overlaps a known segmental duplication.
    """
    for hit in annotations_in_region(var, "segdup", "bed"):
        return True
    return False


def get_conservation_info(var):
    """
    Returns a boolean indicating whether or not the
    variant overlaps a conserved region as defined
    by the 29-way mammalian conservation study.
    http://www.nature.com/nature/journal/v478/n7370/full/nature10530.html

    Data file provenance:
    http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/ \
    29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz

    # Script to convert for gemini:
    gemini/annotation_provenance/make-29way-conservation.sh
    """
    for hit in annotations_in_region(var, "conserved", "bed"):
        return True
    return False


def get_recomb_info(var):
    """
    Returns the mean recombination rate at the site.
    """
    count = 0
    tot_rate = 0.0
    for hit in annotations_in_region(var, "recomb", "bed"):
        if hit.contig not in ['chrY']:
        # recomb rate file is in bedgraph format.
        # pysam will store the rate in the "name" field
            count += 1
            tot_rate += float(hit.name)

    return float(tot_rate) / float(count) if count > 0 else None


def _get_first_vcf_hit(hit_iter):
    if hit_iter is not None:
        hits = list(hit_iter)
        if len(hits) > 0:
            return hits[0]


def _get_vcf_info_attrs(hit):
    info_map = {}
    for info in hit.info.split(";"):
        if info.find("=") > 0:
            (key, value) = info.split("=", 1)
            info_map[key] = value
    return info_map


def get_gms(var):
    """Return Genome Mappability Scores for multiple technologies.
    """
    techs = ["illumina", "solid", "iontorrent"]
    GmsTechs = collections.namedtuple("GmsTechs", techs)
    hit = _get_first_vcf_hit(
        annotations_in_region(var, "gms", "vcf", "grch37"))
    attr_map = _get_vcf_info_attrs(hit) if hit is not None else {}
    return apply(GmsTechs,
                 [attr_map.get("GMS_{0}".format(x), None) for x in techs])


def get_grc(var):
    """Return GRC patched genome regions.
    """
    regions = set()
    for hit in annotations_in_region(var, "grc", "bed", "grch37"):
        regions.add(hit.name)
    return ",".join(sorted(list(regions))) if len(regions) > 0 else None

def get_cse(var):
    """Return if a variant is in a CSE: Context-specific error region.
    """
    for hit in annotations_in_region(var, "cse", "bed", "grch37"):
        return True
    return False

def get_encode_tfbs(var):
    """
    Returns a comma-separated list of transcription factors that were
    observed to bind DNA in this region.  Each hit in the list is constructed
    as TF_CELLCOUNT, where:
      TF is the transcription factor name
      CELLCOUNT is the number of cells tested that had nonzero signals

    NOTE: the annotation file is in BED format, but pysam doesn't
    tolerate BED files with more than 12 fields, so we just use the base
    tuple parser and grab the name column (4th column)
    """
    tfbs = []
    for hit in annotations_in_region(var, "encode_tfbs", "tuple"):
        tfbs.append(hit[3] + "_" + hit[4])
    if len(tfbs) > 0:
        return ','.join(tfbs)
    else:
        return None


def get_encode_dnase_clusters(var):
    """
    If a variant overlaps a DnaseI cluster, return the number of cell types
    that were found to have DnaseI HS at in the given interval, as well
    as a comma-separated list of each cell type:

    Example data:
    chr1	20042385	20042535	4	50.330600	8988t;K562;Osteobl;hTH1
    chr1	20043060	20043210	3	12.450500	Gm12891;T47d;hESCT0
    chr1	20043725	20043875	2	5.948180	Fibrobl;Fibrop
    chr1	20044125	20044275	3	6.437350	HESC;Ips;hTH1
    """
    for hit in annotations_in_region(var, "encode_dnase1", "tuple"):
        return ENCODEDnaseIClusters(hit[3], hit[5])
    return ENCODEDnaseIClusters(None, None)


def get_encode_consensus_segs(var):
    """
    Queries a meta-BEDGRAPH of consensus ENCODE segmentations for 6 cell types:
    gm12878, h1hesc, helas3, hepg2, huvec, k562

    Returns a 6-tuple of the predicted chromatin state of each cell type for the
    region overlapping the variant.

    CTCF: CTCF-enriched element
    E:    Predicted enhancer
    PF:   Predicted promoter flanking region
    R:    Predicted repressed or low-activity region
    TSS:  Predicted promoter region including TSS
    T:    Predicted transcribed region
    WE:   Predicted weak enhancer or open chromatin cis-regulatory element
    """
    for hit in annotations_in_region(var, "encode_consensus_segs", "tuple"):
        return ENCODESegInfo(hit[3], hit[4], hit[5], hit[6], hit[7], hit[8])

    return ENCODESegInfo(None, None, None, None, None, None)


def get_encode_segway_segs(var):
    """
    Queries a meta-BEDGRAPH of SegWay ENCODE segmentations for 6 cell types:
    gm12878, h1hesc, helas3, hepg2, huvec, k562

    Returns a 6-tuple of the predicted chromatin state of each cell type for the
    region overlapping the variant.
    """
    for hit in annotations_in_region(var, "encode_segway_segs", "tuple"):
        return ENCODESegInfo(hit[3], hit[4], hit[5], hit[6], hit[7], hit[8])

    return ENCODESegInfo(None, None, None, None, None, None)


def get_encode_chromhmm_segs(var):
    """
    Queries a meta-BEDGRAPH of SegWay ENCODE segmentations for 6 cell types:
    gm12878, h1hesc, helas3, hepg2, huvec, k562

    Returns a 6-tuple of the predicted chromatin state of each cell type for the
    region overlapping the variant.
    """
    for hit in annotations_in_region(var, "encode_chromhmm_segs", "tuple"):
        return ENCODESegInfo(hit[3], hit[4], hit[5], hit[6], hit[7], hit[8])

    return ENCODESegInfo(None, None, None, None, None, None)


def get_resources( args ):
    """Retrieve list of annotation resources loaded into gemini.
    """
    anno_files = get_anno_files( args )
    return [(n, os.path.basename(anno_files[n])) for n in sorted(anno_files.keys())]
