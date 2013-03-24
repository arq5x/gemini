import zlib
import cPickle
import sqlite3
import numpy

# http://stackoverflow.com/questions/695794/more-efficient-way-to-
# pickle-a-string


def pack_blob(obj):
    return sqlite3.Binary(zdumps(obj))


def unpack_genotype_blob(blob):
    return numpy.array(cPickle.loads(zlib.decompress(blob)))


def zdumps(obj):
    return zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 9)


def zloads(obj):
    return cPickle.loads(zlib.decompress(obj))
