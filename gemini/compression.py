import cPickle
import sqlite3
import numpy
import collections

from gemini_utils import OrderedDict

# http://stackoverflow.com/questions/695794/more-efficient-way-to-
# pickle-a-string

def pack_blob(obj):
    return sqlite3.Binary(zdumps(obj))

def unpack_genotype_blob(blob):
    return cPickle.loads(blob)

def unpack_ordereddict_blob(blob):
    blob_val = cPickle.loads(blob)
    if blob_val is not None:
        return OrderedDict(blob_val)
    return None

def zdumps(obj):
    return cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL)

def zloads(obj):
    return cPickle.loads(obj)
