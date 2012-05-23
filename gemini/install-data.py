#!/usr/bin/env python
import sys
import os
import subprocess
import gemini

"""
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
'ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz.tbi']

def install_annotation_files(anno_root_dir):
    """Download required annotation files.
    """
    
    # create the full gemini data path based on
    # the root dir the user provided
    anno_dir = anno_root_dir + "/gemini/data"

    # update the gemini configuration file
    # with the data directory so that we can 
    # retrieve annotation files at runtime.
    gemini_installation_path = gemini.__path__[0]
    gemini_config_file = os.path.join(gemini_installation_path, 'data/gemini.conf')
    gemini_conf = open(gemini_config_file, 'a')
    gemini_conf.write(anno_dir + '\n')
    gemini_conf.close()

    if not os.path.exists(anno_root_dir + "/gemini"):
        os.mkdir(anno_root_dir + "/gemini")
        os.mkdir(anno_root_dir + "/gemini/data")
    # download and install each of the annotation files
    for dl in anno_files:
        url = "http://people.virginia.edu/~arq5x/files/gemini/annotations/{fname}".format(fname=dl)
        _download_to_dir(url, anno_dir)
    
    
def _download_to_dir(url, dirname):
    """
    Grab an annotation file and place in /usr/share/gemini/data
    """
    print "* downloading " + url + " to " + dirname + "\n"
    stub = os.path.basename(url)
    # download data file to pwd
    cmd = ["curl", "-OL", url]
    subprocess.check_call(cmd)
    # move to system directory (/usr/share/gemini/data) and remove from pwd
    dest = os.path.join(dirname, stub)
    os.rename(stub, dest)

if __name__ == "__main__":
    anno_dir = sys.argv[1]
    if os.path.isdir(anno_dir):
        install_annotation_files(anno_dir)
    else:
        sys.exit(anno_dir + " is not a valid directory.")