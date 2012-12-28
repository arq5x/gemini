#!/usr/bin/env python
import sqlite3
import collections

def map_samples_to_indicies(c):
    """Return a dict mapping samples names (key)
       to sample indices in the numpy genotype arrays (value).
    """
    sample_to_idx = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = str(row['name'])
        idx = row['sample_id'] - 1
        sample_to_idx[name] = idx
    return sample_to_idx


def map_indicies_to_samples(c):
    """Return a dict mapping samples indices in the 
       numpy arrays (key) to sample names.
    """
    idx_to_sample = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = str(row['name'])
        idx = row['sample_id'] - 1
        idx_to_sample[idx] = name
    return idx_to_sample


def get_col_names_and_indices(sqlite_description, ignore_gt_cols = False):
    """Return a list of column namanes and a list of the row indicies.
       Optionally exclude gt_* columns.
    """
    col_indices = []
    col_names = []
    for idx, col_tup in enumerate(sqlite_description):
        # e.g., each col in sqlite desc is a tuple like:
        # ('variant_id', None, None, None, None, None, None)
        col_name = col_tup[0]
        if ((not ignore_gt_cols) or \
           (ignore_gt_cols and not col_name.startswith('gt'))):
            col_indices.append(idx)
            col_names.append(col_name)
    return col_names, col_indices


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

            