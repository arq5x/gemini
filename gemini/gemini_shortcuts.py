#!/usr/bin/env python
import sqlite3
import numpy as np
import zlib
import cPickle

##########################################################################
# Convenience functions
##########################################################################
def map_samples_to_indicies(c):
    """Return a dict mapping samples names (key)
       to sample indices in the numpy genotype arrays (value).
    """
    sample_to_idx = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = row['name']
        idx = row['sample_id'] - 1
        sample_to_idx[name] = idx
    return sample_to_idx


def map_indicies_to_samples(c):
    """Return a dict mapping samples indices in the 
       numpy arrays (key) to sample names.
    """
    idx_to_sample = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = row['name']
        idx = row['sample_id'] - 1
        idx_to_sample[idx] = name
    return idx_to_sample


def get_col_names_and_indices(sqlite_description, ignore_gt_cols = False):
    """Return a list of column namanes and a list of the row indicies.
       Optionally exclude gt_* columns.
    """
    col_indices = []
    col_names = []
    for idx, col_tup in enumerate(sqlite_description):
        # e.g., each col in sqlite desc is a tuple like:
        # ('variant_id', None, None, None, None, None, None)
        col_name = col_tup[0]
        if ((not ignore_gt_cols) or \
           (ignore_gt_cols and not col_name.startswith('gt'))):
            col_indices.append(idx)
            col_names.append(col_name)
    return col_names, col_indices


##########################################################################
# Shortcuts
##########################################################################

def shortcut_variants(args, c):
    """
    Report all columns in the variant table, except for the
    genotype vectors. 
    """
    query = "SELECT * FROM variants"
    c.execute(query)
    
    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        get_col_names_and_indices(c.description, ignore_gt_cols=True)

    if args.use_header:
        print args.separator.join(col for col in col_names)
    for row in c:
        print args.separator.join(str(row[i]) if row[i] is not None else "." \
                                              for i in non_gt_idxs )


def shortcut_genotypes(args, c):
    """For each variant, report each sample's genotype
       on a separate line.
    """
    # grab a mapping of an index in the gts array to the sample name 
    idx_to_sample = map_indicies_to_samples(c)
    
    query = "SELECT  v.chrom, v.start, v.end, \
                     v.ref, v.alt, \
                     v.type, v.sub_type, \
                     v.aaf, v.in_dbsnp, v.gene, \
                     v.gts \
             FROM    variants v"
    c.execute(query)
    
    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        get_col_names_and_indices(c.description, ignore_gt_cols=True)
    col_names.append('sample')
    col_names.append('genotype')
    
    if args.use_header: 
        print args.separator.join(col for col in col_names)
    for row in c:
        gts = np.array(cPickle.loads(zlib.decompress(row['gts'])))
        for idx, gt in enumerate(gts):
            # xrange(len(row)-1) to avoid printing v.gts
            print args.separator.join(str(row[i]) for i in xrange(len(row)-1)),
            print args.separator.join([idx_to_sample[idx], gt])


def shortcut_samples(args, c):
    """
    Report all of the information about the samples in the DB
    """
    query = "SELECT * FROM samples"
    c.execute(query)
    
    (col_names, col_idxs) = get_col_names_and_indices(c.description)
    if args.use_header: 
        print args.separator.join(col_names)
    for row in c:
        print args.separator.join(str(row[i]) if row[i] is not None else "." \
                                              for i in xrange(len(row)) )


def shortcut_snpcounts(args, c):
    """
    Report the count of each type of SNP.
    """
    query = "SELECT ref, alt, count(1) \
           FROM   variants \
           WHERE  type = \'snp\' \
           GROUP BY ref, alt"
    # get the ref and alt alleles for all snps.
    c.execute(query)

    if args.use_header: 
        print args.separator.join(['type', 'count'])
    for row in c:
        print args.separator.join([str(row['ref']) + "->" + str(row['alt']), \
                                   str(row['count(1)'])])


def shortcut_tstv(args, c):
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
    print "transitions" + args.separator + \
          "transversions" + args.separator + "ts/tv"
    print str(ts) + args.separator + \
          str(tv) + args.separator + \
          str(float(ts)/float(tv))


def shortcut_tstv_coding(args, c):
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
    print "transitions" + args.separator + \
          "transversions" + args.separator + "ts/tv"
    print str(ts) + args.separator + \
          str(tv) + args.separator + \
          str(float(ts)/float(tv))


def shortcut_tstv_noncoding(args, c):
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
    print "transitions" + args.separator + \
          "transversions" + args.separator + "ts/tv"
    print str(ts) + args.separator + \
          str(tv) + args.separator + \
          str(float(ts)/float(tv))


def shortcut_sfs(args, c):
    """
    Report the site frequency spectrum
    """
    query = "SELECT round(aaf," + str(args.precision) + "), count(1) \
             FROM (select aaf from variants group by variant_id) \
             GROUP BY round(aaf," + str(args.precision) + ")"
             
    c.execute(query)
    if args.use_header:
        print args.separator.join(['aaf', 'count'])
    for row in c:
        print args.separator.join([str(row[0]), str(row[1])])


def shortcut_mds(c):
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