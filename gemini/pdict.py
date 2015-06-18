"""
This is a ordered dictionary-like object that can handle repeated keys.
It is needed since a query like:

    select v.variant_id, i.variant_id FROM variants v, variant_impacts i where v.variant_id=i.variant_id

will have 'variant_id' twice.

We use it in GeminiQuery to keep the fields in the desired order.
"""
import itertools as it
import sys

def to_json(obj):
    return dict(obj.items())

class PDict(object):
    r"""
    >>> p = PDict([('a', '2'), ('b', 3)])
    >>> p
    {'a': '2', 'b': 3}
    >>> p['a'] = 44
    >>> p
    {'a': 44, 'b': 3}

    >>> p = PDict([('a', '2'), ('b', 3), ('a', 'again'), ('b', 'b-again')])
    >>> p
    {'a': '2', 'b': 3, 'a': 'again', 'b': 'b-again'}
    >>> str(p)
    '2\t3\tagain\tb-again'
    >>> p.keys()
    ['a', 'b', 'a', 'b']
    >>> other = p.copy()
    >>> p['b'] = 555
    >>> p
    {'a': '2', 'b': 555, 'a': 'again', 'b': 'b-again'}

    >>> other
    {'a': '2', 'b': 3, 'a': 'again', 'b': 'b-again'}
    >>> other['a']
    '2'

    >>> import json
    >>> json.dumps(p, default=to_json)
    '{"a": "again", "b": "b-again"}'

    >>> json.dumps(["asdf"], default=to_json)
    '["asdf"]'

    >>> d = PDict()
    >>> d['asdf'] = 1
    >>> d['asff'] = 2

    >>> del d['asdf']
    >>> d
    {'asff': 2}


    """
    __slots__ = ('_keys', '_vals')

    def __init__(self, items=None):
        if items is None:
            self._keys = []
            self._vals = []
        else:
            if isinstance(items, dict):
                self._keys = items.keys()[:]
                self._vals = items.values()[:]
            elif isinstance(items, PDict):
                self._keys = items._keys[:]
                self._vals = items._vals[:]
            elif hasattr(items, "__iter__"):
                self._keys = [i[0] for i in items]
                self._vals = [i[1] for i in items]
            else:
                raise NotImplemented

    def __repr__(self):
        return "{%s}" % ", ".join("%r: %r" % p for p in it.izip(self._keys, self._vals))

    def __str__(self):
        return "\t".join(map(str, self._vals))

    def add(self, key, val):
        self._keys.append(key)
        self._vals.append(val)

    def get(self, key):
        try:
            idx = self._keys.index(key)
        except ValueError:
            return None
        return self._vals[idx]

    def __getitem__(self, key):
        try:
            idx = self._keys.index(key)
        except ValueError:
            raise KeyError(key)
        return self._vals[idx]

    def __setitem__(self, key, val):
        try:
            idx = self._keys.index(key)
        except ValueError:
            self._keys.append(key)
            self._vals.append(val)
            return
        self._vals[idx] = val

    def keys(self):
        return self._keys

    def values(self):
        return self._vals

    def copy(self):
        return self.__class__(self)

    def __iter__(self):
        for k in self._keys:
            yield k

    def __delitem__(self, key):
        idx = self._keys.index(key)
        self._keys = self._keys[:idx] + self._keys[idx + 1:]
        self._vals = self._vals[:idx] + self._vals[idx + 1:]

    def items(self):
        return it.izip(self._keys, self._vals)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
