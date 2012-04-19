#!/usr/bin/env python
import sqlite3
import os
import sys
import numpy as np
import cPickle
import zlib
import collections

import gemini_utils as util

def get_tstv(c, args):
    """
    Report the transition / transversion ratio.
    """
    ts_cmd = "SELECT count(1) \
           FROM  variants \
           WHERE type = \'snp\' \
           AND   sub_type = \'ts\'"
    tv_cmd = "SELECT count(1) \
          FROM  variants v \
          WHERE type = \'snp\' \
          AND   sub_type = \'tv\'"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)
    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "ts" + '\t' + \
          "tv" + '\t' + "ts/tv"
    print str(ts) + '\t' + \
          str(tv) + '\t' + \
          str(float(ts)/float(tv))


def get_tstv_coding(c, args):
    """
    Report the transition / transversion ratio in coding regions.
    """
    ts_cmd = "SELECT count(1) \
           FROM variants v \
           WHERE v.type = \'snp\' \
           AND v.sub_type = \'ts\' \
           AND v.is_coding = 1"
    tv_cmd = "SELECT count(1) \
          FROM variants v \
          WHERE v.type = \'snp\' \
          AND v.sub_type = \'tv\' \
          AND v.is_coding = 1"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)
    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "ts" + '\t' + \
          "tv" + '\t' + "ts/tv"
    print str(ts) + '\t' + \
          str(tv) + '\t' + \
          str(float(ts)/float(tv))


def get_tstv_noncoding(c, args):
    """
    Report the transition / transversion ratio in coding regions.
    """
    ts_cmd = "SELECT count(1) \
           FROM variants v \
           WHERE v.type = \'snp\' \
           AND v.sub_type = \'ts\' \
           AND v.is_coding = 0"
    tv_cmd = "SELECT count(1) \
          FROM variants v \
          WHERE v.type = \'snp\' \
          AND v.sub_type = \'tv\' \
          AND v.is_coding = 0"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)

    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "ts" + '\t' + \
          "tv" + '\t' + "ts/tv"
    print str(ts) + '\t' + \
          str(tv) + '\t' + \
          str(float(ts)/float(tv))


def get_snpcounts(c, args):
    """
    Report the count of each type of SNP.
    """
    query = "SELECT ref, alt, count(1) \
             FROM   variants \
             WHERE  type = \'snp\' \
             GROUP BY ref, alt"
    
    # get the ref and alt alleles for all snps.
    c.execute(query)
    print '\t'.join(['type', 'count'])
    for row in c:
        print '\t'.join([str(row['ref']) + "->" + str(row['alt']), \
                         str(row['count(1)'])])

def get_sfs(c, args):
    """
    Report the site frequency spectrum
    """
    query = "SELECT round(aaf," + str(args.precision) + "), count(1) \
             FROM (select aaf from variants group by variant_id) \
             GROUP BY round(aaf," + str(args.precision) + ")"
             
    c.execute(query)
    if args.use_header:
        print '\t'.join(['aaf', 'count'])
    for row in c:
        print '\t'.join([str(row[0]), str(row[1])])


def shortcut_mds(c, args):
    """
    Compute the pairwise genetic distance between each sample. 
    """
    idx_to_sample = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        idx_to_sample[int(row['sample_id']) - 1] = row['name']

    query = "SELECT DISTINCT v.variant_id, v.gt_types\
    FROM variants v\
    WHERE v.type = 'snp'"
    c.execute(query)

    # keep a list of numeric genotype values
    # for each sample
    genotypes = defaultdict(list)
    for row in c:
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        for idx, type in enumerate(gt_types):
            genotypes[idx_to_sample[idx]].append(type)

            mds = defaultdict(float)
            deno = defaultdict(float)
            # convert the genotype list for each sample
            # to a numpy array for performance.
            # masks stores an array of T/F indicating which genotypes are
            # known (True, [0,1,2]) and unknown (False [-1]). 
            masks = {}
            for sample in genotypes:
                x = np.array(genotypes[sample])
                genotypes[sample] = x
                masks[sample] = \
                np.ma.masked_where(genotypes[sample]>=0, genotypes[sample]).mask
                # compute the euclidean distance for each s1/s2 combination
                # using numpy's vectorized sum() and square() operations.
                # we use the mask arrays to identify the indices of known genotypes
                # for each sample.  by doing a bitwise AND of the mask arrays for the
                # two samples, we have a mask array of variants where __both__ samples
                # were called.
                for s1 in genotypes:
                    for s2 in genotypes:
                        pair = (s1,s2)
                        # which variants have known genotypes for both samples?
                        both_mask = masks[s1] & masks[s2]
                        gt1 = genotypes[s1]
                        gt2 = genotypes[s2]
                        eucl_dist = float(np.sum(np.square((gt1-gt2)[both_mask]))) \
                        / \
                        float(np.sum(both_mask))
                        mds[pair] = eucl_dist
                        deno[pair] = np.sum(both_mask)

                        for pair in mds:
                            print "\t".join([str(pair), str(mds[pair]/deno[pair])])


def get_variants_by_sample(c, args):
    """
    Report the number of variants observed for each sample
    """
    idx_to_sample = util.map_indicies_to_samples(c)
    
    query = "SELECT * FROM variants"
    c.execute(query)
    
    # TO DO: counting variants per sample would be much faster
    #        if we stored a sample_variants table.
    sample_counts = collections.defaultdict(int)
    for row in c:
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        for idx, gt_type in enumerate(gt_types):
            if gt_type > 0: 
                sample_counts[idx_to_sample[idx]] += 1
    
    # report the samples in DESC'ing order of the number of variants observed.
    print '\t'.join(['sample', 'num_variants'])
    for sample in sorted(sample_counts, key=sample_counts.get, reverse=True):
        print sample + "\t" + str(sample_counts[sample])


def get_gtcounts_by_sample(c, args):
    """
    Report the count of each genotype class
    observed for each sample.
    """
    idx_to_sample = util.map_indicies_to_samples(c)
    
    # light struct for tracking the count of each genotype type foreach sample
    class GenoCounts(object):
        """struct for counting the occ. of each gt type."""
        def __init__(self):
            self.hom_ref = 0
            self.het = 0
            self.hom_alt = 0
            self.unknown = 0

    sample_counts = collections.defaultdict(GenoCounts)
    query = "SELECT * FROM variants"
    c.execute(query)
    # count the number of each genotype type obs. for each sample.
    for row in c:
        gt_types  = np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
        for idx, gt_type in enumerate(gt_types):
            sample = idx_to_sample[idx]
            if gt_type == 0: 
                sample_counts[sample].hom_ref += 1
            elif gt_type == 1: 
                sample_counts[sample].het += 1
            elif gt_type == 2: 
                sample_counts[sample].hom_alt += 1
            elif gt_type == -1: 
                sample_counts[sample].unknown += 1

    # report.
    print '\t'.join(['sample', 'num_hom_ref', 'num_het', 
                     'num_hom_alt', 'num_unknown'])
    for sample in sorted(sample_counts, key=sample_counts.get, reverse=True):
        print "\t".join(str(s) for s in [sample, 
                                         sample_counts[sample].hom_ref,
                                         sample_counts[sample].het,
                                         sample_counts[sample].hom_alt,
                                         sample_counts[sample].hom_alt])


def stats(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        if args.tstv:
            get_tstv(c, args)
        elif args.tstv_coding:
            get_tstv_coding(c, args)
        elif args.tstv_noncoding:
            get_tstv_noncoding(c, args)
        elif args.snp_counts:
            get_snpcounts(c, args)
        elif args.sfs:
            get_sfs(c, args)
        elif args.variants_by_sample:
            get_variants_by_sample(c, args)
        elif args.genotypes_by_sample:
            get_gtcounts_by_sample(c, args)
        elif args.mds:
            get_mds(c, args)

