#!/usr/bin/env python
from __future__ import absolute_import
from . import compression as Z

from . import database
from . import gemini_utils as util
from .GeminiQuery import GeminiQuery
import sqlalchemy as sql


def get_variants(conn, metadata, args):
    """
    Report all columns in the variant table, except for the
    genotype vectors.
    """
    query = "SELECT * FROM variants \
             ORDER BY chrom, start"
    res = conn.execute(query)

    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        util.get_col_names_and_indices(metadata.tables["variants"], ignore_gt_cols=True)

    if args.use_header:
        print(args.separator.join(col for col in col_names))
    for row in res:
        print(args.separator.join('.' if (row[i] is None) else
                                  row[i].encode('utf-8') if type(row[i]) is
                                  unicode else str(row[i]) for i in
                                  non_gt_idxs))


def get_genotypes(conn, metadata, args):
    """For each variant, report each sample's genotype
       on a separate line.
    """
    idx_to_sample = util.map_indices_to_samples(metadata)

    query = "SELECT  v.chrom, v.start, v.end, \
                     v.ref, v.alt, \
                     v.type, v.sub_type, \
                     v.aaf, v.in_dbsnp, v.gene, \
                     v.gts \
             FROM    variants v \
             ORDER BY chrom, start"
    res = conn.execute(sql.text(query))

    # build a list of all the column indices that are NOT
    # gt_* columns.  These will be the columns reported
    (col_names, non_gt_idxs) = \
        util.get_col_names_and_indices(metadata.tables["variants"], ignore_gt_cols=True)
    col_names.append('sample')
    col_names.append('genotype')

    if args.use_header:
        print(args.separator.join(col for col in col_names))

    unpack = Z.unpack_genotype_blob
    import zlib


    for row in res:
        try:
            gts = unpack(row['gts'])
        except zlib.error:
            unpack = Z.snappy_unpack_blob
            gts = unpack(row['gts'])

        for idx, gt in enumerate(gts):
            # range(len(row)-1) to avoid printing v.gts
            a = args.separator.join(str(row[i]) for i in range(len(row)-1))
            b = args.separator.join([idx_to_sample[idx], gt])
            print(args.separator.join((a, b)))

def get_samples(conn, metadata, args):
    """
    Report all of the information about the samples in the DB
    """
    query = "SELECT * FROM samples"
    res = conn.execute(query)

    (col_names, col_idxs) = util.get_col_names_and_indices(metadata.tables["samples"])
    if args.use_header:
        print(args.separator.join(col_names))
    for row in res:
        print(args.separator.join(str(row[i]) if row[i] is not None else "." \
                                              for i in range(len(row)) ))


def tfam(args):
    """
    Report the information about the samples in the DB in TFAM format:
    http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    """

    query = ("select family_id, name, paternal_id, maternal_id, "
             "sex, phenotype from samples")
    gq = GeminiQuery(args.db)
    gq.run(query)
    for row in gq:
        print(" ".join(map(str, [row['family_id'], row['name'], row['paternal_id'],
                        row['maternal_id'], row['sex'], row['phenotype']])))


def dump(parser, args):

    conn, metadata = database.get_session_metadata(args.db)

    if args.variants:
        get_variants(conn, metadata, args)
    elif args.genotypes:
        get_genotypes(conn, metadata, args)
    elif args.samples:
        get_samples(conn, metadata, args)
    elif args.tfam:
        tfam(args)
