#!/usr/bin/env python
import argparse
import sys
import os
import shutil
import subprocess
import time

from gemini.config import read_gemini_config, write_gemini_config

"""Install annotation data and update Gemini configuration with location.

The recommended Gemini install location is /usr/local/share/gemini.

    This installation script was inspired by a helpful suggestion from
    Brad Chapman, and is based on his code at:
    https://github.com/chapmanb/bcbb/blob/master/nextgen/tests/test_automated_analysis.py
"""

anno_files = \
['dbsnp.138.vcf.gz',
'clinvar_20140303.vcf.gz',
'29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed.gz',
'hg19.CpG.bed.gz',
'hg19.pfam.ucscgenes.bed.gz',
'hg19.gerp.elements.bed.gz',
'hg19.cytoband.bed.gz',
'hg19.dgv.bed.gz',
'hg19.gwas.bed.gz',
'hg19.rmsk.bed.gz',
'hg19.segdup.bed.gz',
'ESP6500SI.all.snps_indels.vcf.gz',
'ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.2012Oct12.vcf.gz',
'genetic_map_HapMapII_GRCh37.gz',
'wgEncodeRegTfbsClusteredV2.cell_count.20130213.bed.gz',
'encode.6celltypes.consensus.bedg.gz',
'stam.125cells.dnaseI.hg19.bed.gz',
'GRCh37-gms-mappability.vcf.gz',
'GRC_patch_regions.bed.gz',
'kegg_pathways_ensembl66',
'kegg_pathways_ensembl67',
'kegg_pathways_ensembl68',
'kegg_pathways_ensembl69',
'kegg_pathways_ensembl70',
'kegg_pathways_ensembl71',
'hprd_interaction_graph',
'cse-hiseq-8_4-2013-02-20.bed.gz',
'hg19.vista.enhancers.20131108.bed.gz',
'hg19.cosmic.v67.20131024.gz',
'detailed_gene_table_v75',
'summary_gene_table_v75',
'cancer_gene_census.20140120.tsv'
]
extra_anno_files = {"gerp_bp": "hg19.gerp.bw", "cadd_score": "whole_genome_SNVs.tsv.compressed.gz"}

toadd_anno_files = []

anno_versions = {
    "GRCh37-gms-mappability.vcf.gz": 2,
    "dbsnp.138.vcf.gz": 2,
    "clinvar_20140303.vcf.gz": 3,
    "hg19.rmsk.bed.gz": 2,
    "detailed_gene_table_v75": 2,
    "summary_gene_table_v75": 2}

def install_annotation_files(anno_root_dir, dl_files=False, extra=None):
    """Download required annotation files.
    """
    # create the full gemini data path based on
    # the root dir the user provided
    if anno_root_dir.endswith(("gemini/data", "gemini/data/", "gemini_data")):
        anno_dir = anno_root_dir
    elif anno_root_dir.endswith(("gemini", "gemini/")):
        anno_dir = os.path.join(anno_root_dir, "data")
    else:
        anno_dir = os.path.join(anno_root_dir, "gemini", "data")

    cur_config = read_gemini_config(allow_missing=True)
    cur_config["annotation_dir"] = os.path.abspath(anno_dir)
    cur_config["versions"] = anno_versions
    write_gemini_config(cur_config)

    if dl_files:
        if not os.path.exists(anno_dir):
            os.makedirs(anno_dir)
        if not os.path.isdir(anno_dir):
            sys.exit(anno_dir + " is not a valid directory.")
        _check_dependencies()
        to_dl = anno_files[:]
        if extra:
            to_dl += [extra_anno_files[x] for x in extra]
        _download_anno_files("https://s3.amazonaws.com/gemini-annotations",
                             to_dl, anno_dir, cur_config)
    #_download_anno_files("https://s3.amazonaws.com/chapmanb/gemini",
    #                     toadd_anno_files, cur_config)

def _check_dependencies():
    """Ensure required tools for download are present.
    """
    print "Checking required dependencies..."
    for cmd, url in [("curl", "http://curl.haxx.se/")]:
        try:
            retcode = subprocess.call([cmd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except OSError:
            retcode = 127
        if retcode == 127:
            raise OSError("gemini requires %s (%s)" % (cmd, url))
        else:
            print " %s found" % cmd

def _download_anno_files(base_url, file_names, anno_dir, cur_config):
    """Download and install each of the annotation files
    """
    for orig in file_names:
        if orig.endswith(".gz"):
            dls = [orig, "%s.tbi" % orig]
        else:
            dls = [orig]
        for dl in dls:
            url = "{base_url}/{fname}".format(fname=dl, base_url=base_url)
            _download_to_dir(url, anno_dir, anno_versions.get(orig, 1),
                             cur_config.get("versions", {}).get(orig, 1))

def _download_to_dir(url, dirname, version, cur_version):
    """
    Grab an annotation file and place in /usr/share/gemini/data
    """
    print "* downloading " + url + " to " + dirname + "\n"
    stub = os.path.basename(url)
    dest = os.path.join(dirname, stub)
    if not os.path.exists(dest) or version > cur_version:
        # download data file to staging directory instead of current
        # direction which may not have space
        dl_dir = os.path.join(dirname, "tmpdownload")
        if not os.path.exists(dl_dir):
            os.makedirs(dl_dir)
        orig_dir = os.getcwd()
        os.chdir(dl_dir)
        # download file, allowing retries for network errors
        max_retries = 2
        retries = 0
        while 1:
            cmd = ["curl", "-C", "-", "-OL", url]
            retcode = subprocess.call(cmd)
            if retcode == 0:
                break
            else:
                print "Curl failed with non-zero exit code %s. Retrying" % retcode
                if retries >= max_retries:
                    raise
                time.sleep(10)
                retries += 1
        with open(stub) as in_handle:
            line1 = in_handle.readline()
            line2 = in_handle.readline()
            if "?xml" in line1 and "Error" in line2 and "AccessDenied" in line2:
                raise ValueError("Could not download annotation file. Permission denied error: %s" % url)
        # move to system directory (/usr/share/gemini/data) and remove from pwd
        shutil.move(stub, dest)
        os.chdir(orig_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("anno_dir", help="Directory to write annotation files.")
    parser.add_argument("--nodata", dest="dl_files", default=True, action="store_false")
    parser.add_argument("--extra", help="Add additional non-standard genome annotations to include",
                        action="append", default=[])
    args = parser.parse_args()
    install_annotation_files(args.anno_dir, args.dl_files, args.extra)
