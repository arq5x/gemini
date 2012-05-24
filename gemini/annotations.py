#!/usr/bin/env python
import pysam
import os
import sys
import collections

# determine where the user installed the gemini annotation suite.
gemini_installation_path = os.path.split(__file__)[0]
gemini_config_file = os.path.join(gemini_installation_path, 'data/gemini.conf')
gemini_conf = open(gemini_config_file, 'r')
config_lines = gemini_conf.readlines()
if len(config_lines) > 0:
    anno_dirname = config_lines[0].rstrip()
else:
    sys.exit('Cannot determine where gemini annotation files are located.  Exiting.')
    
# dictionary of anno_type -> open Tabix file handles
annos = {}
def load_annos():
    """
    Populate a dictionary of Tabixfile handles for
    each annotation file.  Other modules can then
    access a given handle and fetch data from it
    as follows:
    
    dbsnp_handle = annotations.annos['dbsnp']
    hits = dbsnp_handle.fetch(chrom, start, end)
    """    
    anno_files   = {
                    'cytoband'  : os.path.join(anno_dirname, 'hg19.cytoband.bed.gz'),
                    'dbsnp'     : os.path.join(anno_dirname, 'dbsnp.135.vcf.gz'),
                    'gwas'      : os.path.join(anno_dirname, 'hg19.gwas.bed.gz'),
                    'rmsk'      : os.path.join(anno_dirname, 'hg19.rmsk.bed.gz'),
                    'segdup'    : os.path.join(anno_dirname, 'hg19.segdup.bed.gz'),
                    'conserved' : os.path.join(anno_dirname, '29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed.gz'),
                    'cpg_island': os.path.join(anno_dirname, 'hg19.CpG.bed.gz'),
                    'dgv'       : os.path.join(anno_dirname, 'hg19.dgv.bed.gz'),
                    'esp'       : os.path.join(anno_dirname, 'ESP5400.all.snps.vcf.gz'),
                    '1000g'     : os.path.join(anno_dirname, 'ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz'),
                    'recomb'    : os.path.join(anno_dirname, 'genetic_map_HapMapII_GRCh37.gz')
                   }

    for anno in anno_files:
        annos[anno] = pysam.Tabixfile(anno_files[anno])

DbSnpInfo = collections.namedtuple("DbSnpInfo", "rs_ids in_omim clin_sig")

def get_cpg_island_info(var):
    """
    Returns a boolean indicating whether or not the
    variant overlaps a CpG island 
    """
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    for hit in annos['cpg_island'].fetch(chrom, var.start, var.end, parser=pysam.asBed()):
        return True
    return False


def get_cyto_info(var):
    """
    Returns a comma-separated list of the chromosomal
    cytobands that a variant overlaps.
    """
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    cyto_band = ''
    for hit in annos['cytoband'].fetch(chrom, var.start, var.end, 
                                       parser=pysam.asBed()): 
        if len(cyto_band) > 0:
            cyto_band += "," + chrom + hit.name
        else: 
            cyto_band += chrom + hit.name
    return cyto_band if len(cyto_band) > 0 else None


def get_dbsnp_info(var):
    """
    Returns a suite of annotations from dbSNP
    """
    chrom = var.CHROM if not var.CHROM.startswith("chr") else var.CHROM[3:]
    rs_ids  = []
    clin_sigs = []
    in_omim = 0
    for hit in annos['dbsnp'].fetch(chrom, var.start, var.end, parser=pysam.asVCF()):
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


def get_esp_info(var):
    """
    Returns a suite of annotations from the ESP project
    """
    ESPInfo = collections.namedtuple("ESPInfo", "found aaf_EA aaf_AA aaf_ALL exome_chip")
     
    chrom = var.CHROM if not var.CHROM.startswith("chr") else var.CHROM[3:]
    aaf_EA = aaf_AA = aaf_ALL = None
    maf = fetched = con = []
    exome_chip = False
    found = False
    info_map = {}
    if chrom not in ['Y']:
        for hit in annos['esp'].fetch(chrom, var.start, var.end, parser=pysam.asVCF()):
            fetched.append(hit)
            # We need a single ESP entry for a variant
            if fetched != None and len(fetched) == 1 and hit.alt == var.ALT[0] and hit.ref == var.REF:
                found = True    
                # loads each VCF INFO key/value pair into a DICT
                for info in hit.info.split(";"):      
                    if info.find("=") > 0:
                    # splits on first occurence of '='   
                    # useful to handle valuerror: too many values to unpack (e.g (a,b) = split(",", (a,b,c,d)) for cases like
                    #SA=http://www.ncbi.nlm.nih.gov/sites/varvu?gene=4524&amp%3Brs=1801131|http://omim.org/entry/607093#0004  
                        (key, value) = info.split("=", 1)
                        info_map[key] = value
                # get the % minor allele frequencies      
                if info_map.get('MAF') is not None:
                    lines = info_map['MAF'].split(",")
                    # divide by 100 because ESP reports allele frequencies as percentages.
                    aaf_EA = float(lines[0]) / 100.0
                    aaf_AA = float(lines[0]) / 100.0
                    aaf_ALL = float(lines[0]) / 100.0
                    
                #Is the SNP on an human exome chip?
                if info_map.get('EXOME_CHIP') is not None and info_map['EXOME_CHIP'] == "no":
                    exome_chip = 0
                elif info_map.get('EXOME_CHIP') is not None and info_map['EXOME_CHIP'] == "yes":
                    exome_chip = 1        
    return ESPInfo(found, aaf_EA, aaf_AA, aaf_ALL, exome_chip)  
        
def get_1000G_info(var):
    """
    Returns a suite of annotations from the 1000 Genomes project
    """
    ThousandGInfo = collections.namedtuple("ThousandGInfo", 
                                           "found aaf_ALL aaf_AMR aaf_ASN aaf_AFR aaf_EUR")

    chrom = var.CHROM if not var.CHROM.startswith("chr") else var.CHROM[3:]
    fetched = []
    info_map = {}
    found = False
    for hit in annos['1000g'].fetch(chrom, var.start, var.end, parser=pysam.asVCF()):
        fetched.append(hit)
        # We need a single 1000G entry for a variant
        if fetched != None and len(fetched) == 1 and hit.alt == var.ALT[0] and hit.ref == var.REF:
            # loads each VCF INFO key/value pair into a DICT
            found = True
            for info in hit.info.split(";"):      
                if info.find("=") > 0:
                    (key, value) = info.split("=", 1)
                    info_map[key] = value

    return ThousandGInfo(found, info_map.get('AF'), info_map.get('AMR_AF'), 
                         info_map.get('ASN_AF'), info_map.get('AFR_AF'), 
                         info_map.get('EUR_AF'))
            
def get_rmsk_info(var):
    """
    Returns a comma-separated list of annotated repeats
    that overlap a variant.  Derived from the UCSC rmsk track 
    """
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    rmsk_hits = []
    for hit in annos['rmsk'].fetch(chrom, var.start, var.end, parser=pysam.asBed()):
        rmsk_hits.append(hit.name)
    return ",".join(rmsk_hits) if len(rmsk_hits) > 0 else None


def get_segdup_info(var):
    """
    Returns a boolean indicating whether or not the
    variant overlaps a known segmental duplication. 
    """
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    for hit in annos['segdup'].fetch(chrom, var.start, var.end, parser=pysam.asBed()):
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
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    for hit in annos['conserved'].fetch(chrom, var.start, var.end, parser=pysam.asBed()):
        return True
    return False

def get_recomb_info(var):
    """
    Returns the mean recombination rate at the site.
    """
    chrom = var.CHROM if var.CHROM.startswith("chr") else "chr" + var.CHROM
    count = 0
    tot_rate = 0.0
    if chrom not in ['chrY']:
        # recomb rate file is in bedgraph format.
        # pysam will store the rate in the "name" field
        for hit in annos['recomb'].fetch(chrom, var.start, var.end, 
                                           parser=pysam.asBed()): 
            count += 1
            tot_rate += float(hit.name)

    return float(tot_rate) / float(count) if count > 0 else None
