#!/usr/bin/env python
import sqlite3
import os
import sys
import numpy as np
import cPickle
import zlib
import collections
from copy import copy

import gemini_utils as util
from gemini_constants import *

class Site(object):
    def __init__(self, row):
        self.chrom     = row['chrom']
        self.start     = row['start']
        self.end       = row['end']
        self.ref       = row['ref']
        self.alt       = row['alt']
        self.gene      = row['gene']
        self.num_hets  = row['num_het']
        self.exon      = row['exon']
        self.aaf       = row['aaf']
        self.in_dbsnp  = row['in_dbsnp']
        self.impact    = row['impact']
        self.is_lof    = row['is_lof']
        # specific to each sample
        self.phased    = None
        self.gt        = None
    
    def __repr__(self):
        return ','.join([self.chrom, str(self.start), str(self.end),
                         self.ref, self.alt, self.gt,
                         self.impact, self.exon, str(self.aaf), 
                         str(self.in_dbsnp)])
        
    def __eq__(self, other):
        return self.start == other.start


def get_compound_hets(c, args):
    """
    Report candidate compound heterozygous mutations.
    """
    # build a mapping of the numpy array index to the appropriate sample name
    # e.g. 0 == 109400005
    #     37 == 147800025
    idx_to_sample = util.map_indicies_to_samples(c)

    comp_hets = collections.defaultdict(lambda: collections.defaultdict(list))

    query = "SELECT * FROM variants \
             WHERE is_coding = 1" # is_exonic - what about splice?
    c.execute(query)

    
    # step 1. collect all candidate heterozygptes for all
    # genes and samples.  the list will be refined in step 2.
    for row in c:
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        gt_phases = np.array(cPickle.loads(zlib.decompress(row['gt_phases'])))
        gt_bases  = np.array(cPickle.loads(zlib.decompress(row['gts'])))
        
        site = Site(row)
        
        # filter putative sites that the user doesn't care about
        if site.num_hets > 1 and not args.allow_other_hets: 
            continue
        if not site.is_lof and args.only_lof:
            continue
        
        # track each sample that is heteroyzgous at this site.
        for idx, gt_type in enumerate(gt_types):
            if gt_type == GT_HET:
                sample = idx_to_sample[idx]
                # (testing)
                # sample = "NA19002"
                sample_site = copy(site)
                sample_site.phased = gt_phases[idx]
                
                # require phased genotypes
                if not sample_site.phased:
                    continue
                sample_site.gt = gt_bases[idx]
                # add the site to the list of candidates
                # for this sample/gene
                comp_hets[sample][site.gene].append(sample_site)
    
    
    # header
    print "sample\tgene\thet1\thet2"
    # step 2.  now, cull the list of candidate heterozygotes for each
    # gene/sample to those het pairs where the alternate alleles
    # were inherited on opposite haplotypes.
    for sample in comp_hets:
        for gene in comp_hets[sample]:
            for site1 in comp_hets[sample][gene]:
                for site2 in comp_hets[sample][gene]:
                    if site1 == site2:
                        continue
                    
                    # expand the genotypes for this sample
                    # at each site into it's composite
                    # alleles.  e.g. A|G -> ['A', 'G']
                    alleles_site1 = site1.gt.split('|')
                    alleles_site2 = site2.gt.split('|')
                    
                    # return the haplotype on which the alternate
                    # allele was observed for this sample at each
                    # candidate het. site.
                    # e.g., if ALT=G and alleles_site1=['A', 'G']
                    # then alt_hap_1 = 1.  if ALT=A, then alt_hap_1 = 0
                    alt_hap_1 = alleles_site1.index(site1.alt)
                    alt_hap_2 = alleles_site2.index(site2.alt)
                    
                    # it is only a true compound heterozygote iff
                    # the alternates are on opposite haplotypes.
                    if alt_hap_1 != alt_hap_2:
                        print "\t".join([sample, gene, str(site1), str(site2)])


def run(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        # run the compound het caller
        get_compound_hets(c, args)

