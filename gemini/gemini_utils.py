#!/usr/bin/env python
import sqlite3

def map_samples_to_indicies(c):
    """Return a dict mapping samples names (key)
       to sample indices in the numpy genotype arrays (value).
    """
    sample_to_idx = {}
    c.execute("select sample_id, name from samples")
    for row in c:
        name = row['name']
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
        name = row['name']
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