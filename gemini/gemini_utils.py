#!/usr/bin/env python
from __future__ import absolute_import

import numpy as np
import collections
from collections import defaultdict
try:
    from itertools import tee, ifilterfalse
    PY3 = False
except ImportError:
    from itertools import tee, filterfalse as ifilterfalse
    basestring = str
    PY3 = True

import locale
ENC = locale.getpreferredencoding()

def to_str(s, enc=ENC):
    if hasattr(s, "decode"):
        return s.decode(enc)
    if isinstance(s, np.str_):
        return str(s)
    return s

from .gemini_subjects import Subject
import sqlalchemy as sql

def get_gt_cols(metadata):
    try:
        return [c.name for c in metadata.tables["variants"].columns if
                c.name.startswith("gt") and
                c.type.__class__.__name__.upper() == "BLOB"]
    except KeyError:
        # if there's no variants table, there are no gt_cols.
        return []

def map_indices_to_samples(metadata):
    """Return a dict mapping samples indices in the
       numpy arrays (key) to sample names.
    """
    samples = list(metadata.tables["samples"].select().execute())
    d = {s['sample_id'] - 1: s['name'] for s in samples}
    return [d[i] for i in range(len(d))]

def get_col_names_and_indices(tbl, ignore_gt_cols=False):
    """Return a list of column names and a list of the row indices.
       Optionally exclude gt_* columns.
    """
    inames = [(i, c) for i, c in enumerate(tbl.columns)]
    if ignore_gt_cols:
        inames = [x for x in inames if not (x[1].startswith("gt") and
                  x[1].type.__class__.__name__.upper() == "BLOB")]

    return [x[1].name for x in inames], [x[0] for x in inames]

# http://code.activestate.com/recipes/576694/
class OrderedSet(collections.MutableSet):

    def __init__(self, iterable=None):
        self.end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

try:
    from thread import get_ident as _get_ident
except ImportError:
    try:
        from threading import get_ident as _get_ident
    except ImportError:
        from dummy_thread import get_ident as _get_ident

try:
    from _abcoll import KeysView, ValuesView, ItemsView
except ImportError:
    pass


def itersubclasses(cls, _seen=None):
    """
    snagged from:  http://code.activestate.com/recipes/576949/
    itersubclasses(cls)

    Generator over all subclasses of a given class, in depth first order.

    >>> list(itersubclasses(int)) == [bool]
    True
    >>> class A(object): pass
    >>> class B(A): pass
    >>> class C(A): pass
    >>> class D(B,C): pass
    >>> class E(D): pass
    >>>
    >>> for cls in itersubclasses(A):
    ...     print(cls.__name__)
    B
    D
    E
    C
    >>> # get ALL (new-style) classes currently defined
    >>> [cls.__name__ for cls in itersubclasses(object)] #doctest: +ELLIPSIS
    ['type', ...'tuple', ...]
    """

    if not isinstance(cls, type):
        raise TypeError('itersubclasses must be called with '
                        'new-style classes, not %.100r' % cls)
    if _seen is None:
        _seen = set()
    try:
        subs = cls.__subclasses__()
    except TypeError:  # fails only when cls is type
        subs = cls.__subclasses__(cls)
    for sub in subs:
        if sub not in _seen:
            _seen.add(sub)
            yield sub
            for sub in itersubclasses(sub, _seen):
                yield sub

def partition(pred, iterable):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = tee(iterable)
    return list(ifilterfalse(pred, t1)), list(filter(pred, t2))

def quote_string(item):
    """ if the item is a string, put quotes around it else leave it """
    if isinstance(item, basestring):
        item = "\"" + item + "\""
    return item

def partition_by_fn(seq, key_fn=lambda x: x, val_fn=lambda x: x):
    """
    partition a sequence into a dictionary with keys key_fn(x) and
    list of values val_fn(x)
    """
    d = defaultdict(list)
    for x in seq:
        d[key_fn(x)].append(val_fn(x))
    return d
