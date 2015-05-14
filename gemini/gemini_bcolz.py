"""
This is an example module for using genotype-query engines for gemini.
An engine only has to have a query() function that returns a list of
variant ids that meet a gemini genotype query.

Any engine can be added in a post-hoc fashion provided a module is
provided with the single filter():

    # take a gemini --gt-filter string and return a list of variant_ids that meet
    # that filter. `obj` is the thing returned by load().
    # user_dict contains things like HET, UNKNOWN, etc. used in the eval.
    filter(db_path, obj, gt_filter, user_dict)

See below for an implementation using bcolz.
It is using carray rather than ctable because with ctable, we'd be limited to
2908 samples on a ext3 file-system. with carrays, we can have up to 31998
samples on ext3. These limits are not an issue for ext4.

"""

import os
import sys
sys.setrecursionlimit(8192)
import sqlite3
import time
import re
import shutil

import numpy as np
import bcolz
import numexpr as ne
bcolz.blosc_set_nthreads(2)
ne.set_num_threads(2)

import compression
from gemini_utils import get_gt_cols
decomp = compression.unpack_genotype_blob

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

def create(db, cols=None):
    if cols is None:
        cols = [x[0] for x in gt_cols_types if x[0] != 'gts']
        print >>sys.stderr, (
                "indexing all columns execpt 'gts'; to index that column, "
                "run gemini bcolz_index %s --cols gts" % db)

    conn = sqlite3.connect(db)
    cur = conn.cursor()
    gt_cols = [x for x in get_gt_cols(cur) if x in cols]
    samples = get_samples(cur)
    bcpath = get_bcolz_dir(db)

    mkdir(bcpath)

    nv = get_n_variants(cur)

    sys.stderr.write("loading %i variants for %i samples into bcolz\n"
                     % (nv, len(samples)))

    if nv == 0 or len(samples) == 0:
        return

    carrays = {}
    tmps = {}
    try:
        for gtc in gt_cols:
            carrays[gtc] = []
            tmps[gtc] = []

            dt = dict(gt_cols_types)[gtc]
            for s in samples:
                mkdir("%s/%s" % (bcpath, s))
                carrays[gtc].append(bcolz.carray(np.empty(0, dtype=dt),
                                    expectedlen=nv,
                                    rootdir="%s/%s/%s" % (bcpath, s, gtc),
                                    chunklen=16384*8,
                                    mode="w"))
                tmps[gtc].append([])

        t0 = time.time()
        step = 200000
        del gtc

        empty = [-1] * len(samples)
        for i, row in enumerate(cur.execute("select %s from variants" % ", ".join(gt_cols))):
            for j, gt_col in enumerate(gt_cols):
                vals = decomp(row[j])
                if vals is None: # empty gt_phred_ll
                    vals = empty
                for isamp, sample in enumerate(samples):
                    tmps[gt_col][isamp].append(vals[isamp])
                    if (i > 0 and i % step == 0) or i == nv - 1:
                        carrays[gt_col][isamp].append(tmps[gt_col][isamp])
                        tmps[gt_col][isamp] = []
                        carrays[gt_col][isamp].flush()

            if i % step == 0 and i > 0:
                print >>sys.stderr, "at %.1fM (%.0f rows / second)" % (i / 1000000., i / float(time.time() - t0))

        t = float(time.time() - t0)
        print >>sys.stderr, "loaded %d variants at %.1f / second" % (len(carrays[gt_col][0]), nv / t)
    except:
        # on error, we remove the dirs so we can't have weird problems.
        for k, li in carrays.items():
            for i, ca in enumerate(li):
                if i < 5:
                    print >>sys.stderr, "removing:", ca.rootdir
                if i == 5:
                    print >>sys.stderr, "not reporting further removals for %s" % k
                ca.flush()
                shutil.rmtree(ca.rootdir)
        raise


# TODO: since we call this from query, we can improve speed by only loading
# samples that appear in the query with an optional query=None arg to load.
def load(db):

    t0 = time.time()
    conn = sqlite3.connect(db)
    cur = conn.cursor()

    gt_cols = get_gt_cols(cur)
    samples = get_samples(cur)
    bcpath = get_bcolz_dir(db)

    carrays = {}
    for gtc in gt_cols:
        carrays[gtc] = []
        for s in samples:
            path = "%s/%s/%s" % (bcpath, s, gtc)
            if os.path.exists(path):
                carrays[gtc].append(bcolz.open(path, mode="r"))
    if os.environ.get("GEMINI_DEBUG") == "TRUE":
        print >>sys.stderr, "it took %.2f seconds to load arrays" \
            % (time.time() - t0)
    return carrays


def filter(db, query, user_dict):
    # these should be translated to a bunch or or/and statements within gemini
    # so they are supported, but must be translated before getting here.
    if "any(" in query or "all(" in query or \
       ("sum(" in query and not query.startswith("sum(") and query.count("sum(") == 1):
        return None
    user_dict['where'] = np.where

    carrays = load(db)
    if query.startswith("not "):
        # "~" is not to numexpr.
        query = "~" + query[4:]
    sum_cmp = False
    if query.startswith("sum("):
        assert query[-1].isdigit()
        query, sum_cmp = query[4:].rsplit(")", 1)
        query = "(%s) %s" % (query, sum_cmp)

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
    if os.environ.get('GEMINI_DEBUG') == 'TRUE':
        print >>sys.stderr, query

    # loop through and create a cache of "$gt__$sample"
    for gt_col in carrays:
        # if not gt_col in query: continue
        for i, sample_array in enumerate(carrays[gt_col]):
            sample = samples[i]
            # if not sample in query: continue
            user_dict["%s__%s" % (gt_col, sample)] = sample_array

    # had to special-case count. it won't quite be as efficient
    if "|count|" in query:
        tokens = query[2:-2].split("|count|")
        icmp = tokens[-1]
        # a list of carrays, so not in memory.
        res = [bcolz.eval(tok, user_dict=user_dict) for tok in tokens[:-1]]
        # in memory after this, but just a single axis array.
        res = np.sum(res, axis=0)
        res = ne.evaluate('res%s' % icmp)
    else:
        res = bcolz.eval(query, user_dict=user_dict)

    if res.shape[0] == 1 and len(res.shape) > 1:
        res = res[0]
    variant_ids, = np.where(res)
    #variant_ids = np.array(list(bcolz.eval(query, user_dict=user_dict,
    #    vm="numexpr").wheretrue()))
    # variant ids are 1-based.
    if len(variant_ids) > 0:
        return 1 + variant_ids
    else:
        return []

if __name__ == "__main__":

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

