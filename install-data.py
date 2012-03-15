#!/usr/bin/env python
import sys
import os
import subprocess

def install_annotation_files():
    """Download required annotation files.
    """
    anno_files = \
    ['dbsnp.135.vcf.gz', 'dbsnp.135.vcf.gz.tbi',
    'hg19.29way.bed.gz', 'hg19.29way.bed.gz.tbi',
    'hg19.CpG.bed.gz', 'hg19.CpG.bed.gz.tbi',
    'hg19.cytoband.bed.gz', 'hg19.cytoband.bed.gz.tbi',
    'hg19.dgv.bed.gz', 'hg19.dgv.bed.gz.tbi',
    'hg19.gwas.bed.gz', 'hg19.gwas.bed.gz.tbi',
    'hg19.rmsk.bed.gz', 'hg19.rmsk.bed.gz.tbi',
    'hg19.segdup.bed.gz', 'hg19.segdup.bed.gz.tbi']

    dirname = "/usr/share/pop/data/"
    if not os.path.exists(dirname):
        os.mkdir("/usr/share/pop")
        os.mkdir("/usr/share/pop/data")
    # download and install each of the annotation files
    for dl in anno_files:
        url = "http://people.virginia.edu/~arq5x/files/pop/annotations/{fname}".format(fname=dl)
        _download_to_dir(url, dirname)

def _download_to_dir(url, dirname):
    """
    Grab an annotation file and place in /usr/share/pop/data
    """
    print "* downloading " + url + " to " + dirname + "\n"
    stub = os.path.basename(url)
    # download data file to pwd
    cl = ["curl", "-OL", url]
    subprocess.check_call(cl)
    # move to system directory (/usr/share/pop/data) and remove from pwd
    dest = os.path.join(dirname, stub)
    os.rename(stub, dest)

if __name__ == "__main__":
    install_annotation_files()
