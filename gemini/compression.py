import zlib
import cPickle

try:
    from cyordereddict import OrderedDict
except:
    from collections import OrderedDict


def pack_blob(obj):
    return zdumps(obj)

def unpack_genotype_blob(blob):
    return cPickle.loads(zlib.decompress(blob))

def unpack_ordereddict_blob(blob):
    blob_val = cPickle.loads(zlib.decompress(blob))
    if blob_val is not None:
        return OrderedDict(blob_val)
    return None

def zdumps(obj):
    return zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 9)

def zloads(obj):
    return cPickle.loads(zlib.decompress(obj))
