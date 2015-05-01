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
    # user_dict contains things like HET, UNKNOWN, etc. used in the eval.
    query(db_path, obj, gt_filter, user_dict)

See below for an implementation using bcolz.
It is using carray rather than ctable because with ctable, we'd be limited to
2908 samples on a ext3 file-system. with carrays, we can have up to 31998
samples on ext3. These limits are not an issue for ext4.

"""

import os
import sys
import sqlite3
import time
import re

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

def create(db, cols=[x[0] for x in gt_cols_types]):
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    gt_cols = [x for x in get_gt_cols(cur) if x in cols]
    samples = get_samples(cur)
    bcpath = get_bcolz_dir(db)

    mkdir(bcpath)

    nv = get_n_variants(cur)

    print >>sys.stderr, "loading %i variants for %i samples" % (nv, len(samples))

    carrays = {}
    tmps = {}
    for gtc in gt_cols:
        carrays[gtc] = []
        tmps[gtc] = []

        dt = dict(gt_cols_types)[gtc]
        for s in samples:
            mkdir("%s/%s" % (bcpath, s))
            carrays[gtc].append(bcolz.carray(np.empty(0, dtype=dt),
                expectedlen=nv, rootdir="%s/%s/%s" % (bcpath, s, gtc),
                mode="w"))
            tmps[gtc].append([])


    t0 = time.time()

    empty = [-1] * len(samples)
    for i, row in enumerate(cur.execute("select %s from variants" % ", ".join(gt_cols))):
        for j, gt_col in enumerate(gt_cols):
            vals = decomp(row[j])
            if vals is None: # empty gt_phred_ll
                vals = empty
            for isamp, sample in enumerate(samples):
                tmps[gt_col][isamp].append(vals[isamp])
                if i % 20000 == 0 or i == nv - 1:
                    carrays[gt_col][isamp].append(tmps[gt_col][isamp])
                    tmps[gt_col][isamp] = []

        if i % 20000 == 0:
            print >>sys.stderr, "at %i" % i

    t = float(time.time() - t0)
    print >>sys.stderr, "loaded %d variants at %.1f / second" % (len(carrays[gt_col][isamp]), nv / t)

def load(db):
    conn = sqlite3.connect(db)
    cur = conn.cursor()

    gt_cols = get_gt_cols(cur)
    samples = get_samples(cur)
    bcpath = get_bcolz_dir(db)

    carrays = {}
    for gtc in gt_cols:
        carrays[gtc] = []
        for s in samples:
            carrays[gtc].append(
                    bcolz.open("%s/%s/%s" % (bcpath, s, gtc), mode="r"))
    return carrays

def query(db, carrays, query, user_dict):
    if "any(" in query or "all(" in query:
        return None
    if carrays is None:
        carrays = load(db)
    query = query.replace(".", "__")
    query = " & ".join("(%s)" % token for token in query.split(" and "))
    query = " | ".join("(%s)" % token for token in query.split(" or "))

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    samples = get_samples(cur)

    # convert gt_col[index] to gt_col__sample_name
    patt = "(%s)\[(\d+)\]" % "|".join(carrays.keys())

    def subfn(x):
        """Turn gt_types[1] into gt_types__sample"""
        field, idx = x.groups()
        return "%s__%s" % (field, samples[int(idx)])

    query = re.sub(patt, subfn, query)
    print >>sys.stderr, query

    # loop through and create a cache of "$gt__$sample"
    for gt_col in carrays:
        # if not gt_col in query: continue
        for i, sample_array in enumerate(carrays[gt_col]):
            sample = samples[i]
            # if not sample in query: continue
            user_dict["%s__%s" % (gt_col, sample)] = sample_array

    print [x for x in user_dict.keys() if not x.startswith("gt")]
    print user_dict['sample_info']
    variant_ids, = np.where(bcolz.eval(query, user_dict=user_dict, vm="numexpr"))
    # variant ids are 1-based.
    if len(variant_ids) > 0:
        return 1 + variant_ids
    else:
        return []

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

    print query(db, carrays, q, user_dict=dict(HET=1, HOM_REF=0, HOM_ALT=3,
        UNKNOWN=2))
    print "compare to:", ("""gemini query -q "select variant_id, gts.1719PC0016 from variants" """
                          """ --gt-filter "%s" %s""" % (q, db))

