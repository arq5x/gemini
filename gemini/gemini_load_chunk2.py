#!/usr/bin/env python
from __future__ import print_function
from effects import SnpEff

# native Python imports
import os.path
import time
import sys
import sqlite3
import itertools as it

import toml  # toml.py

# third-party imports
import cyvcf2 as vcf
import blosc
blosc.set_nthreads(1)
blosc.set_blocksize(8192)

import zlib
import cPickle

def opack_blob(obj, _none=buffer(zlib.compress(cPickle.dumps(None, cPickle.HIGHEST_PROTOCOL)))):
    if obj is None: return _none
    return buffer(zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 1))

def pack_blob(obj):
    if obj is None: return ''
    return buffer(blosc.compress(obj.tostring(), obj.dtype.itemsize, clevel=5, shuffle=True))
    #return buffer(blosc.pack_array(obj))

def is_number(op, field):
    return field.endswith("_float") or op in ("mean", "median", "min", "max")

def is_integer(op, field):
    return field.endswith("_int")

def is_flag(op, infos):
    if not op in infos: return False
    return infos[op].type == "Flag"

def name_type(op, field, infos):
    """
    >>> name_type("sum", "t_float", {})
    ('t', 'REAL')
    >>> name_type("aaa", "tt_int", {})
    ('tt', 'INTEGER')
    >>> name_type("aaa", "ingt_flag", {})
    ('ingt', 'BOOL')
    >>> name_type("aaa", "hello", {})
    ('hello', 'TEXT')
    >>> name_type("flag", "LCR", {})
    ('LCR', 'BOOL')
    """
    if field.endswith("_float"):
        return field[:-6], "REAL"
    if field.endswith("_int"):
        return field[:-4], "INTEGER"
    if field.endswith("flag"):
        return field[:-5], "BOOL"
    if is_flag(field, infos):
        return field, "BOOL"
    if op == "flag":
        return field, "BOOL"
    if op in ("min", "mean", "max"):
        return field, "REAL"

    return field, "TEXT"

class VCFLoader(object):
    def __init__(self, path, toml_path, ped=None, db_path="test.db"):

        # this will find the fields we are expecting during load.
        self.vcf_reader = vcf.VCF(path)
        self.load_config(toml_path)
        self.ped = ped
        if db_path is None:
            db_path = (path[:-6 if path.endswith('.vcf.gz') else -3]) + "db"
        self.db_path = db_path
        self.db = sqlite3.connect(self.db_path)
        self.db.isolation_level = None
        self.cursor = self.db.cursor()
        self.cursor.execute('PRAGMA synchronous = OFF')
        self.cursor.execute('PRAGMA journal_mode=MEMORY')

        self.create_tables()

        self.extended_fields = ["chrom", "start", "end", "variant_id", "vid", "ref",
                                "alt", "qual", "filter", "type", "sub_type"] \
                               + self.fields + ["gts", "gt_types", "gt_phases",
                                       "gt_depths", "gt_ref_depths",
                                       "gt_alt_depths", "gt_quals",
                                       #"gt_copy_numbers",
                                       "gt_phred_ll_homref", "gt_phred_ll_het", "gt_phred_ll_homalt"]
    def create_tables(self):
        self._create_variants_table()
        self._create_sample_table()

    def _create_sample_table(self):
        samples = self.vcf_reader.samples
        self.sample_to_idx = {s: i for i, s in enumerate(samples, start=1)}
        create = """\
CREATE TABLE samples (
    sample_id INTEGER,
    family_id text default NULL,
    name text default NULL,
    paternal_id text default NULL,
    maternal_id TEXT default NULL,
    sex TEXT default NULL,
    phenotype TEXT default NULL%s
)"""

        q = "INSERT INTO samples VALUES (%s)" % ",".join(["?"] * 7)

        if self.ped is None:
            self.cursor.execute(create)
            self.cursor.executemany(q, [(self.sample_to_idx[s], 0, s, 0, 0, -9, -9) for s in samples])
            return

        fh = open(self.ped)
        lheader = next(fh)
        header = ["sample_id", "family_id", "paternal_id", "maternal_id", "sex", "phenotype"]
        if lheader[0] != "#":
            fh = it.chain([lheader], fh)
        else:
            header = header + lheader[1:].strip().split()[6:]

        create = create % (",\n    " + ",\n    ".join(("%s TEXT" % h) for h in header[6:]))
        self.cursor.execute(create)

        inserts = []
        for toks in (l.split("\t") if l.count("\t") > 1 else l.split() for l in fh):
            # set name as sample_id
            toks.insert(2, toks[1])
            # switch family_id and sample_id since they are in opposite order in # ped
            toks[:2] = toks[:2][::-1]
            try:
                toks[0] = self.sample_to_idx[toks[0]]
            except KeyError:
                # present in ped but not in samples
                continue
            inserts.append(toks)

        assert len(set(len(x) for x in inserts)) == 1, ("weird ped file. missing tabs and spaces?")
        q = "INSERT INTO samples VALUES (%s)" % ",".join(["?"] * len(inserts[0]))
        self.cursor.executemany(q, sorted(inserts))

    def _create_variants_table(self):
        tmpl = """\
CREATE table variants (
    chrom TEXT,
    start INTEGER,
    end INTEGER,
    variant_id INTEGER,
    vid TEXT,
    ref TEXT,
    alt TEXT,
    qual REAL,
    filter TEXT,
    type TEXT,
    sub_type TEXT,
    gene TEXT,
    impact_severity TEXT,
    is_exonic BOOL,
    is_coding BOOL,
    is_lof BOOL,
    biotype TEXT,
    sift_pred TEXT,
    sift_score REAL,
    polyphen_pred TEXT,
    polyphen_score REAL,
    %s,
    gts BLOB,
    gt_types BLOB,
    gt_phases BLOB,
    gt_depths BLOB,
    gt_ref_depths BLOB,
    gt_alt_depths BLOB,
    gt_quals BLOB,
    --gt_copy_numbers BLOB,
    gt_phred_ll_homref BLOB,
    gt_phred_ll_het BLOB,
    gt_phred_ll_homalt BLOB,
    PRIMARY KEY(variant_id ASC)
    )"""
        cmd = tmpl % (",\n    ".join("%s %s" % (field, ftype) for field, ftype in
                      it.izip(self.fields, self.types)))
        self.cursor.execute(cmd)

    def load_config(self, toml_path):
        cfg = toml.loads(open(toml_path).read())
        fields = []
        ops = []
        for f in cfg['annotation']:
            fields.extend(f.get('names', f.get('fields')))
            ops.extend(f['ops'])
        assert len(ops) == len(fields)
        self.config = cfg
        self.o_fields = fields
        self.ops = ops

        field_types = [name_type(o, f, {}) for f, o in it.izip(self.o_fields, self.ops)]
        self.fields, self.types = map(list, zip(*field_types))

    def insert(self, variants):
        if len(variants) == 0:
            return
        t0 = time.time()
        cur = self.cursor
        query = "INSERT INTO variants VALUES(%s)" % ",".join(["?"] * len(variants[0]))
        try:
            cur.executemany(query, variants)
        except sqlite3.ProgrammingError:
            raise
            for var in variants:
                cur.execute(query, var)
        print("insert time:", time.time() - t0)

    def load(self, buffer_size=10000):
        tl = time.time()
        load_buffer = []
        ots = zip(self.types, self.o_fields)
        tgt, igt, bt = 0, 0, 0
        # TODO: pull the columns from the table and use a dict. To hard to track
        # 140 columns.
        for i, v in enumerate(self.vcf_reader, start=1):
            t0 = time.time()
            vals = [v.INFO_get(f, '' if t == 'TEXT' else False if t == 'BOOL' else  None) for t, f in ots]
            igt += time.time() - t0

            t0 = time.time()
            gts = [
                opack_blob(v.gt_bases),
                pack_blob(v.gt_types),
                pack_blob(v.gt_phases),
                pack_blob(v.gt_depths),
                pack_blob(v.gt_ref_depths),
                pack_blob(v.gt_alt_depths),
                pack_blob(v.gt_quals),
                #pack_blob(np.array(v.gt_copy_numbers, np.float32)),
                pack_blob(v.gt_phred_ll_homref),
                pack_blob(v.gt_phred_ll_het),
                pack_blob(v.gt_phred_ll_homalt)]
            tgt += time.time() - t0

            t0 = time.time()
            ann = SnpEff.top_severity(v.INFO_get('ANN').split(","))
            if isinstance(ann, list):
                ann = ann[0]

            basic = [v.CHROM, v.start, v.end, i - 1, v.ID or '', v.REF,
                     ",".join(v.ALT), v.QUAL, v.FILTER, v.var_type,
                     v.var_subtype, ann.gene, ann.impact_severity,
                     ann.exonic, ann.coding, ann.lof, ann.biotype,
                     ann.sift_class, ann.sift_value,
                     ann.polyphen_class, ann.polyphen_value]

            bt += time.time() - t0
            variant = basic + vals + gts
            load_buffer.append(variant)

            if i % buffer_size == 0:
                self.insert(load_buffer)
                load_buffer = []
        self.insert(load_buffer)
        print("gt time", tgt)
        print("info time", igt)
        print("basic time", bt)
        print("total load time:", time.time() - tl)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    v = VCFLoader(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    v.load()
