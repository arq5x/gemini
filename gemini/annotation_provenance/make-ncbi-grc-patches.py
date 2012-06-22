#!/usr/bin/env python
"""Retrieve information on patches and fixes to GRCh37 from GRC website.

Converts information on these regions into a BED file for Gemini import.

http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
"""
import urllib2
from collections import namedtuple
from operator import attrgetter
from contextlib import closing

patch = "p8"
base_url = "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/" \
           "vertebrates_mammals/Homo_sapiens/GRCh37.{0}/".format(patch)
sub_url = "/alt_scaffolds/alt_scaffold_placement.txt"
dirs = ["PATCHES"] + ["ALT_REF_LOCI_%s" % i for i in range(1, 10)]

def main():
    out_file = "GRC_patch_regions.bed"
    regions = []
    for dname in dirs:
        cur_url = base_url + dname + sub_url
        for region in grc_regions_from_url(cur_url):
            regions.append(region)
    regions.sort(key=attrgetter("chrom", "start", "end"))
    with open(out_file, "w") as out_handle:
        for region in regions:
            out_handle.write("{chrom}\t{start}\t{end}\t{name}\n".format(
                    **vars(region)))

def grc_regions_from_url(url):
    GrcRegion = namedtuple("GrcRegion", "chrom,start,end,name")
    with closing(urllib2.urlopen(url)) as in_handle:
        header = in_handle.next()
        for parts in (l.split("\t") for l in in_handle):
            try:
                chrom = int(parts[5])
            except ValueError:
                chrom = parts[5]
            yield GrcRegion(chrom, int(parts[11]) - 1, int(parts[12]),
                            "grc_%s" % ("fix" if parts[2].endswith("PATCH") else "novel"))
        

if __name__ == "__main__":
    main()
