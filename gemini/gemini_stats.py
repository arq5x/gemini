#!/usr/bin/env python
import numpy as np
import collections
from collections import Counter
import compression as Z

import sqlalchemy as sql
import gemini_utils as util
from gemini_constants import *
import GeminiQuery


def get_tstv(conn, metadata, args):
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
    res = conn.execute(sql.text(ts_cmd))
    ts = res.fetchone()[0]
    # get the number of transversions
    res = conn.execute(sql.text(tv_cmd))
    tv = res.fetchone()[0]
    # report the transitions, transversions, and the ts/tv ratio
    print "ts" + '\t' + \
          "tv" + '\t' + "ts/tv"
    print str(ts) + '\t' + \
        str(tv) + '\t' + \
        str(tstv(ts, tv))


def get_tstv_coding(conn, metadata, args):
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
    res = conn.execute(ts_cmd)
    ts = res.fetchone()[0]

    # get the number of transversions
    res = conn.execute(tv_cmd)
    tv = res.fetchone()[0]

    # report the transitions, transversions, and the ts/tv ratio
    print "ts" + '\t' + \
          "tv" + '\t' + "ts/tv"
    print str(ts) + '\t' + \
        str(tv) + '\t' + \
        str(tstv(ts, tv))


def get_tstv_noncoding(conn, metadata, args):
    """
    Report the transition / transversion ratio in non-coding regions.
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
    res = conn.execute(ts_cmd)
    ts = res.fetchone()[0]

    # get the number of transversions
    res = conn.execute(tv_cmd)
    tv = res.fetchone()[0]

    # report the transitions, transversions, and the ts/tv ratio
    print "ts" + '\t' + \
          "tv" + '\t' + "ts/tv"
    print str(ts) + '\t' + \
        str(tv) + '\t' + \
        str(tstv(ts, tv))


def tstv(ts, tv):
    """
    Calculate ts/tv, and avoid division by zero error
    """
    try:
        return round(float(ts) / float(tv), 4)
    except ZeroDivisionError:
        return 0

def get_snpcounts(conn, metadata, args):
    """
    Report the count of each type of SNP.
    """
    query = "SELECT ref, alt, count(1) \
             FROM   variants \
             WHERE  type = \'snp\' \
             GROUP BY ref, alt"

    # get the ref and alt alleles for all snps.
    res = conn.execute(sql.text(query))
    print '\t'.join(['type', 'count'])
    for row in res:
        print '\t'.join([str(row['ref']) + "->" + str(row['alt']),
                         str(row['count(1)'])])


def get_sfs(conn, metadata, args):
    """
    Report the site frequency spectrum
    """
    precision = 3
    query = "SELECT round(aaf," + str(precision) + "), count(1) \
             FROM variants \
             GROUP BY round(aaf," + str(precision) + ")"

    res = conn.execute(sql.text(query))
    print '\t'.join(['aaf', 'count'])
    for row in res:
        print '\t'.join([str(row[0]), str(row[1])])


def get_mds(conn, metadata, args):
    """
    Compute the pairwise genetic distance between each sample.
    """
    idx_to_sample = {}
    res = conn.execute(sql.text("select sample_id, name from samples"))
    for row in res:
        idx_to_sample[int(row['sample_id']) - 1] = row['name']

    query = "SELECT DISTINCT v.variant_id, v.gt_types\
    FROM variants v\
    WHERE v.type = 'snp'"
    res = conn.execute(query)

    # keep a list of numeric genotype values
    # for each sample
    genotypes = collections.defaultdict(list)
    import zlib
    unpack = Z.unpack_genotype_blob

    for row in res:

        try:
            gt_types = unpack(row['gt_types'])
        except zlib.error:
            unpack = Z.snappy_unpack_blob
            gt_types = unpack(row['gt_types'])

        # at this point, gt_types is a numpy array
        # idx:  0 1 2 3 4 5 6 .. #samples
        # type [0 1 2 1 2 0 0 ..         ]
        for idx, gt_type in enumerate(gt_types):
            sample = idx_to_sample[idx]
            genotypes[sample].append(gt_type)

    mds = collections.defaultdict(float)
    # convert the genotype list for each sample
    # to a numpy array for performance.
    # masks stores an array of T/F indicating which genotypes are
    # known (True, [0,1,2]) and unknown (False [-1]).
    masks = {}
    for s in genotypes:
        sample = str(s)
        x = np.array(genotypes[sample])
        genotypes[sample] = x
        masks[sample] = \
            np.ma.masked_where(genotypes[sample] != UNKNOWN,
                               genotypes[sample]).mask

    # compute the euclidean distance for each s1/s2 combination
    # using numpy's vectorized sum() and square() operations.
    # we use the mask arrays to identify the indices of known genotypes
    # for each sample.  by doing a bitwise AND of the mask arrays for the
    # two samples, we have a mask array of variants where __both__ samples
    # were called.
    for sample1 in genotypes:
        for sample2 in genotypes:
            pair = (sample1, sample2)
            # which variants have known genotypes for both samples?
            both_mask = masks[str(sample1)] & masks[str(sample2)]
            genotype1 = genotypes[sample1]
            genotype2 = genotypes[sample2]

            # distance between s1 and s2:
            eucl_dist = float(np.sum(np.square((genotype1 - genotype2)[both_mask]))) \
                / \
                float(np.sum(both_mask))

            mds[pair] = eucl_dist

    # report the pairwise MDS for each sample pair.
    print "sample1\tsample2\tdistance"
    for pair in mds:
        print "\t".join([str(pair[0]), str(pair[1]), str(round(mds[pair], 4))])


def get_variants_by_sample(conn, metadata, args):
    """
    Report the number of variants observed for each sample
    where the sample had a non-ref genotype
    """
    idx_to_sample = util.map_indices_to_samples(metadata)

    # report.
    print '\t'.join(['sample', 'total'])

    query = "SELECT sample_id, \
             (num_het + num_hom_alt) as total \
             FROM sample_genotype_counts"
    res = conn.execute(sql.text(query))
    for row in res:
        sample = idx_to_sample[row['sample_id'] - 1]
        print "\t".join(str(s) for s in [sample,
                                         row['total']])


def get_gtcounts_by_sample(conn, metadata, args):
    """
    Report the count of each genotype class
    observed for each sample.
    """
    idx_to_sample = util.map_indices_to_samples(metadata)

    # report.
    print '\t'.join(['sample', 'num_hom_ref', 'num_het',
                     'num_hom_alt', 'num_unknown', 'total'])

    query = "SELECT *, \
             (num_hom_ref + num_het + num_hom_alt + num_unknown) as total \
             FROM sample_genotype_counts"
    res = conn.execute(query)
    # count the number of each genotype type obs. for each sample.
    for row in res:
        sample = idx_to_sample[row['sample_id'] -1]
        print "\t".join(str(s) for s in [sample,
                                         row['num_hom_ref'],
                                         row['num_het'],
                                         row['num_hom_alt'],
                                         row['num_unknown'],
                                         row['total']])


def summarize_query_by_sample(args):
    gq = GeminiQuery.GeminiQuery(args.db)
    gq.run(args.query, show_variant_samples=True, gt_filter=args.gt_filter)
    total_counts = Counter()
    het_counts = Counter()
    hom_alt_counts = Counter()
    hom_ref_counts = Counter()
    print "\t".join(["sample", "total", "num_het", "num_hom_alt", "num_hom_ref"])
    for row in gq:
        total_counts.update(row["variant_samples"])
        het_counts.update(row["het_samples"])
        hom_alt_counts.update(row["hom_alt_samples"])
        hom_ref_counts.update(row["hom_ref_samples"])
    for key in total_counts.keys():
        count_row = [key, total_counts.get(key, 0), het_counts.get(key, 0),
                     hom_alt_counts.get(key, 0), hom_ref_counts.get(key, 0)]
        print "\t".join(map(str, count_row))


def stats(parser, args):

    import database
    conn, metadata = database.get_session_metadata(args.db)

    if args.tstv:
        get_tstv(conn, metadata, args)
    elif args.tstv_coding:
        get_tstv_coding(conn, metadata, args)
    elif args.tstv_noncoding:
        get_tstv_noncoding(conn, metadata, args)
    elif args.snp_counts:
        get_snpcounts(conn, metadata, args)
    elif args.sfs:
        get_sfs(conn, metadata, args)
    elif args.variants_by_sample:
        get_variants_by_sample(conn, metadata, args)
    elif args.genotypes_by_sample:
        get_gtcounts_by_sample(conn, metadata, args)
    elif args.mds:
        get_mds(conn, metadata, args)
    elif args.query:
        summarize_query_by_sample(args)
