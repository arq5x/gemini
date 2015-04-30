"""
This is an example module for using genotype-query engines for gemini.
An engine only has to have a query() function that returns a list of
variant ids that meet a gemini genotype query.

Any engine can be added in a post-hoc fashion provided a module is
provided with the following functions:

    # create the index given the database for the first time
    create(db_path)
    # load an existing genotype-query engine (this returns any object needed to query).
    # this may not be needed by some engines and can return None in that case.
    load(db_path)
    # take a gemini --gt-filter string and return a list of variant_ids that meet
    # that filter. `obj` is the thing returned by load().
    query(db_path, obj, gt_filter)

See below for an implementation using bcolz.
It is using carray rather than ctable because with ctable, we'd be limited to
2908 samples on a ext3 file-system. with carrays, we can have up to 31998
samples on ext3. These limits are not an issue for ext4.

"""

import sqlite3
import time
import os

import numpy as np
import bcolz

import compression
decomp = compression.unpack_genotype_blob

def get_gt_cols(cur):
    keys = ('col_num', 'col_name', 'col_type', '_', '_', '_')
    gts = []
    for row in cur.execute("pragma table_info(variants)"):
        d = dict(zip(keys, row))
        if d['col_name'][:2] == 'gt' and d['col_type'] == 'blob':
            gts.append(d['col_name'])
    return gts

def get_samples(cur):
    return [x[0] for x in cur.execute("select name from samples")]

def get_n_variants(cur):
    return next(cur.execute("select count(*) from variants"))[0]

def get_bcolz_dir(db):
    return db + ".gts"

gt_cols_types = (
    ('gts', np.object),
    ('gt_types', np.int8),
    ('gt_phases', np.bool),
    ('gt_depths', np.int32),
    ('gt_ref_depths', np.int32),
    ('gt_alt_depths', np.int32),
    ('gt_quals', np.float32),
    ('gt_copy_numbers', np.int32),
    ('gt_phred_ll_homref', np.int32),
    ('gt_phred_ll_het', np.int32),
    ('gt_phred_ll_homalt', np.int32),
)

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass

def create(db):
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    gt_cols = get_gt_cols(cur)
    samples = get_samples(cur)
    bcpath = get_bcolz_dir(db)

    mkdir(bcpath)

    nv = get_n_variants(cur)

    print >>sys.stderr, "loading %i variants for %i samples" % (nv, len(samples))

    carrays = {}
    for name in gt_cols:
        carrays[name] = []
        dt = dict(gt_cols_types)[name]
        for s in samples:
            mkdir("%s/%s" % (bcpath, s))
            carrays[name].append(bcolz.carray(np.empty(0, dtype=dt),
                expectedlen=nv, rootdir="%s/%s/%s" % (bcpath, s, name), mode='w'))

    t0 = time.time()

    for i, row in enumerate(cur.execute("select %s from variants" % ", ".join(gt_cols))):
        for j, gt_col in enumerate(gt_cols):
            vals = decomp(row[j])
            for isamp, sample in enumerate(samples):
                carrays[gt_col][isamp].append(vals[isamp])
                if i % 100 == 0 or i == nv - 1:
                    carrays[gt_col][isamp].flush()

    t = float(time.time() - t0)
    print >>sys.stderr, "loaded %d variants at %.1f / second" % (len(carrays[gt_col][isamp]), nv / t)

def load(db):
    conn = sqlite3.connect(db)
    cur = conn.cursor()

    gt_cols = get_gt_cols(cur)
    samples = get_samples(cur)
    bcpath = get_bcolz_dir(db)

    carrays = {}
    for name in gt_cols:
        carrays[name] = []
        for s in samples:
            carrays[name].append(
                    bcolz.open("%s/%s/%s" % (bcpath, s, name), mode="r"))
    return carrays

def query(db, carrays, query):
    query = query.replace(".", "__")
    query = " & ".join("(%s)" % token for token in query.split(" and "))
    query = " | ".join("(%s)" % token for token in query.split(" or "))

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    samples = get_samples(cur)

    cache = dict(HOM_REF=0, HET=1, UNKNOWN=2, HOM_ALT=3,
                 MISSING=None, UNAFFECTED=1, AFFECTED=2)
    # loop through and create a cache of "$gt__$sample"
    for gt_col in carrays:
        # if not gt_col in query: continue
        for i, sample_array in enumerate(carrays[gt_col]):
            sample = samples[i]
            # if not sample in query: continue
            cache["%s__%s" % (gt_col, sample)] = sample_array

    variant_ids, = np.where(bcolz.eval(query, user_dict=cache, vm="numexpr"))
    # variant ids are 1-based.
    return 1 + variant_ids

if __name__ == "__main__":

    import sys
    db = sys.argv[1]
    #create(sys.argv[1])
    carrays = load(db)
    conn = sqlite3.connect(db)
    if len(sys.argv) > 2:
        q = sys.argv[2]
    else:
        q = "gt_types.1094PC0012 == HET and gt_types.1719PC0016 == HET and gts.1094PC0012 == 'A/C'"

    print query(db, carrays, q)
    print "compare to:", ("""gemini query -q "select variant_id, gts.1719PC0016 from variants" """
                          """ --gt-filter "%s" %s""" % (q, db))

