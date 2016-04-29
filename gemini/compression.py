"""
>>> a = np.arange(20)
>>> (a == snappy_unpack_blob(snappy_pack_blob(a))).all()
True

>>> a = np.arange(20, dtype=np.uint8)
>>> snappy_unpack_blob(snappy_pack_blob(a)).dtype == a.dtype
True

>>> a = np.arange(20, dtype=np.float32)
>>> b = snappy_unpack_blob(snappy_pack_blob(a))
>>> a.dtype == b.dtype == np.float32
True

>>> all(ai == bi for ai, bi in zip(a, b))
True

>>> snappy_unpack_blob(snappy_pack_blob(None)) is None
True

>>> exp = ['A', 'XSDFSDFSD', "A/A,B"]
>>> obs = snappy_unpack_blob(snappy_pack_blob(np.array(exp)))
>>> assert exp == list(obs)

>>> try:
...    b = snappy_pack_blob(a)
... except zlib.error:
...    pass


>>> b = snappy_unpack_blob(snappy_pack_blob(np.array([False, True, False, True])))
>>> b[2] = True
"""
import zlib
import cPickle

try:
    from cyordereddict import OrderedDict
except:
    from collections import OrderedDict


def pack_blob(obj):
    return buffer(zdumps(obj))

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

try:
    import snappy
    import numpy as np
except ImportError:
    pass

# we use the numpy type char as the first item we save to know the dtype when we
# decompress.
lookup = {dt(1).dtype.char: dt for dt in (np.uint8, np.uint32, np.int32,
    np.float32, np.int64, np.float64, np.bool_)}

SEP = '\0'

def snappy_pack_blob(obj, sep=SEP):
    if obj is None: return ''
    c = obj.dtype.char
    if c == 'S': return 'S' + snappy.compress(sep.join(obj))
    return buffer(c + snappy.compress(obj.tobytes()))

def snappy_unpack_blob(blob, sep=SEP):
    if len(blob) == 0: return None
    if blob[0] == 'S':
        return np.array(snappy.decompress(blob[1:]).split(sep))
    dt = lookup[blob[0]]
    arr = np.frombuffer(snappy.decompress(blob[1:]), dtype=dt)
    # hack since arrays arent writable from buffer and we need this for comp_het
    # phasing.
    if blob[0] == '?':
        arr.setflags(write=True)
    return arr

if __name__ == "__main__":
    import doctest
    print(doctest.testmod())
