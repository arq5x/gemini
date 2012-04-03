#!/usr/bin/env python
import pysam
import os
import collections

# pop REQUIRES that the annotation files be installed at /usr/share/pop/data
# This is done at installation time by invoking:
#    % sudo python install-data.py

anno_dirname = '/usr/share/gemini/data/'

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
                    'conserved' : os.path.join(anno_dirname, 'hg19.29way.bed.gz'),
                    'cpg_island': os.path.join(anno_dirname, 'hg19.CpG.bed.gz'),
                    'dgv'       : os.path.join(anno_dirname, 'hg19.dgv.bed.gz')
            }
    for anno in anno_files:
        annos[anno] = pysam.Tabixfile(anno_files[anno])


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
    DbSnpInfo = collections.namedtuple("DbSnpInfo", "rs_ids in_omim clin_sig")

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
