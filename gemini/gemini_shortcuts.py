#!/usr/bin/env python
import sqlite3
import numpy as np
import zlib
import cPickle

def map_samples_to_indicies(c):
    sample_to_idx = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = row['name']
        idx = row['sample_id'] - 1
        sample_to_idx[name] = idx
    return sample_to_idx
    
def map_indicies_to_samples(c):
    idx_to_sample = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = row['name']
        idx = row['sample_id'] - 1
        idx_to_sample[idx] = name
    return idx_to_sample


def shortcut_genotypes(c):
    """For each variant, report each sample's genotype
       on a separate line.
    """
    # grab a mapping of an index in the gts array to the sample name 
    idx_to_sample = map_indicies_to_samples(c)
    
    query = "SELECT  v.chrom, v.start, v.end, \
                     v.ref, v.alt, \
                     v.type, v.sub_type, \
                     v.aaf, v.in_dbsnp, \
                     v.gts \
             FROM    variants v"
    c.execute(query)
    for row in c:
        gts = np.array(cPickle.loads(zlib.decompress(row['gts'])))
        for idx, gt in enumerate(gts):
            print "\t".join(str(row[i]) for i in xrange(9)),
            print "\t".join([idx_to_sample[idx], gt])


def shortcut_ir_candidates(args):
    query = "SELECT v.*, s.name, g.type \
               FROM  variants v, \
                     genotypes g, \
                     samples s \
               WHERE v.variant_id       = g.variant_id \
                AND  g.sample_id       = s.sample_id \
                AND  v.num_hom_alt = 1 \
                AND  v.num_het = 0 \
                AND  g.type   > 0 \
                AND  v.exonic = 1 \
                AND  v.depth >= 200 \
                GROUP BY v.chrom, v.start, v.end, v.impact, v.gene, s.name, g.type" # 1 = het, 2 = hom_alt
    sqlite_cli_call(args, query)


def shortcut_samples(args):
    query = "SELECT * FROM samples"
    sqlite_cli_call(args, query)


def shortcut_snpcounts(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    cmd = "SELECT ref, alt, count(1) \
           FROM   variants \
           WHERE  type = \'snp\' \
           GROUP BY ref, alt"
    # get the ref and alt alleles for all snps.
    c.execute(cmd)

    if args.use_header: args.separator.join(c.description)
    for row in c:
        ref   = str(row[0])
        alt   = str(row[1])
        count = row[2]
        print ref + "->" + alt + args.separator + str(count)


def shortcut_tstv(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
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
    print "transitions" + args.separator + "transversions" + args.separator + "ts/tv"
    print str(ts) + "\t" + str(tv) + "\t" + str(float(ts)/float(tv))


def shortcut_tstv_coding(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    ts_cmd = "SELECT count(1) \
           FROM variants v \
           WHERE v.type = \'snp\' \
           AND v.sub_type = \'ts\' \
           AND v.exonic = 1"
    tv_cmd = "SELECT count(1) \
          FROM variants v \
          WHERE v.type = \'snp\' \
          AND v.sub_type = \'tv\' \
          AND v.exonic = 1"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)
    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "transitions" + args.separator + "transversions" + args.separator + "ts/tv"
    print str(ts) + "\t" + str(tv) + "\t" + str(float(ts)/float(tv))


def shortcut_tstv_noncoding(args):
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    ts_cmd = "SELECT count(1) \
           FROM variants v \
           WHERE v.type = \'snp\' \
           AND v.sub_type = \'ts\' \
           AND v.exonic = 0"
    tv_cmd = "SELECT count(1) \
          FROM variants v \
          WHERE v.type = \'snp\' \
          AND v.sub_type = \'tv\' \
          AND v.exonic = 0"
    # get the number of transitions
    c.execute(ts_cmd)
    ts = c.fetchone()[0]
    # get the number of transversions
    c.execute(tv_cmd)

    tv = c.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "transitions" + args.separator + "transversions" + args.separator + "ts/tv"
    print str(ts) + args.separator + str(tv) + args.separator + str(float(ts)/float(tv))


def shortcut_variants(args):
    query = "SELECT * FROM variants"
    sqlite_cli_call(args, query)


def shortcut_sfs(args):
    query = "SELECT round(aaf," + str(args.precision) + "), count(1) \
             FROM (select aaf from variants group by variant_id) \
             GROUP BY round(aaf," + str(args.precision) +")"
    sqlite_cli_call(args, query)


def shortcut_mds(c):
    """
    Compute the pairwise genetic distance between each sample. 
    """
    idx_to_sample = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        idx_to_sample[int(row['sample_id']) - 1] = row['name']
    
    mds_cmd = "SELECT DISTINCT v.variant_id, v.gt_types\
               FROM variants v\
               WHERE v.type = 'snp'"
    c.execute(mds_cmd)

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