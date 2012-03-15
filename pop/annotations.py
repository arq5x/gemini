#!/usr/bin/env python
import pysam
import os

# pop REQUIRES that the annotation files be installed at /usr/share/pop/data
# This is done at installation time by invoking:
#    % sudo python install-data.py

anno_dirname = '/usr/share/pop/data/'

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
    # dictionary of anno_type -> open Tabix file handles
    annos = {}
    for anno in anno_files:
        annos[anno] = pysam.Tabixfile(anno_files[anno])
    return annos
