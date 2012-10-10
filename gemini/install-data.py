#!/usr/bin/env python
import sys
import os
import shutil
import subprocess

from gemini.config import read_gemini_config, write_gemini_config

"""Install annotation data and update Gemini configuration with location.

The recommended Gemini install location is /usr/local/share/gemini.

    This installation script was inspired by a helpful suggestion from
    Brad Chapman, and is based on his code at:
    https://github.com/chapmanb/bcbb/blob/master/nextgen/tests/test_automated_analysis.py
"""

anno_files = \
['dbsnp.135.vcf.gz', 
'dbsnp.135.vcf.gz.tbi',
'29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed.gz', 
'29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed.gz.tbi',
'hg19.CpG.bed.gz', 
'hg19.CpG.bed.gz.tbi',
'hg19.cytoband.bed.gz', 
'hg19.cytoband.bed.gz.tbi',
'hg19.dgv.bed.gz', 
'hg19.dgv.bed.gz.tbi',
'hg19.gwas.bed.gz', 
'hg19.gwas.bed.gz.tbi',
'hg19.rmsk.bed.gz', 
'hg19.rmsk.bed.gz.tbi',
'hg19.segdup.bed.gz', 
'hg19.segdup.bed.gz.tbi',
'ESP5400.all.snps.vcf.gz', 
'ESP5400.all.snps.vcf.gz.tbi',
'ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz', 
'ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz.tbi',
'genetic_map_HapMapII_GRCh37.gz',
'genetic_map_HapMapII_GRCh37.gz.tbi',
'wgEncodeRegTfbsClusteredV2.cell_count.bed.gz',
'wgEncodeRegTfbsClusteredV2.cell_count.bed.gz.tbi',
'encode.6celltypes.consensus.bedg.gz',
'encode.6celltypes.consensus.bedg.gz.tbi',
'encode.6celltypes.segway.bedg.gz',
'encode.6celltypes.segway.bedg.gz.tbi',
'encode.6celltypes.chromhmm.bedg.gz',
'encode.6celltypes.chromhmm.bedg.gz.tbi',
'GRCh37-gms-mappability.vcf.gz',
'GRCh37-gms-mappability.vcf.gz.tbi',
'GRC_patch_regions.bed.gz',
'GRC_patch_regions.bed.gz.tbi',
'kegg_pathways_ensembl66',
'kegg_pathways_ensembl67',
'kegg_pathways_ensembl68',
'hprd_interaction_graph'
]

def install_annotation_files(anno_root_dir):
    """Download required annotation files.
    """
    # create the full gemini data path based on
    # the root dir the user provided
    if anno_root_dir.endswith(("gemini", "gemini/")):
        anno_dir = os.path.join(anno_root_dir, "data")
    else:
        anno_dir = os.path.join(anno_root_dir, "gemini", "data")
    if not os.path.exists(anno_dir):
        os.makedirs(anno_dir)

    cur_config = read_gemini_config(allow_missing=True)
    cur_config["annotation_dir"] = anno_dir
    write_gemini_config(cur_config)

    # download and install each of the annotation files
    for dl in anno_files:
        url =        "http://people.virginia.edu/~arq5x/files/gemini/annotations/{fname}".format(fname=dl)
        _download_to_dir(url, anno_dir)

 
def _download_to_dir(url, dirname):
    """
    Grab an annotation file and place in /usr/share/gemini/data
    """
    print "* downloading " + url + " to " + dirname + "\n"
    stub = os.path.basename(url)
    dest = os.path.join(dirname, stub)
    if not os.path.exists(dest):
        # download data file to pwd
        cmd = ["curl", "-OL", url]
        subprocess.check_call(cmd)
        # move to system directory (/usr/share/gemini/data) and remove from pwd
        shutil.move(stub, dest)

if __name__ == "__main__":
    anno_dir = sys.argv[1]
    if not os.path.exists(anno_dir):
        os.makedirs(anno_dir)
    if os.path.isdir(anno_dir):
        install_annotation_files(anno_dir)
    else:
        sys.exit(anno_dir + " is not a valid directory.")
