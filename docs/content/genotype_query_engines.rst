#########################
Speeding genotype queries
#########################

Design of Gemini
================

Gemini stores all genotypes information in a single column. E.g. `gt_depths`
is a single column in the sqlite database that contains information on all 
samples in a compressed and serialized array. So, in order to do a query involving
`--gt-filter` `gemini` must iterate over each row, decompress and de-serialize the
array and then evaluate the genotype filter. Even if the filter involves only a
single sample, we must deserialize the entire array. This design works quite well
and we have improved performance greatly, but for complex queries, it is quite slow.

We aim to provide the means to *optionally** provide an additional index that can
be (again) *optionally* used to perform the genotype filtering.
These are external to gemini in that they do not change the behavior of `gemini`
when used without the engine, but they could create additional tables in the
gemini database.

Design of Genotype Query Engines
================================

Genotype Query Engines can be plugged in to `gemini`. They must be
exposed with a single function:

    query(db_path, gt_filter, user_dict)

where `db_path` is the path to the gemini sqlite database, `gt_filter` is
the genotype query string. user_dict will be pre-filled with things like
user_dict contains things like HET, UNKNOWN, etc. used in gemini.

The `query` function must return a list of variant_ids that meet the specified
filter. If it can not perform the query, it must return `None`.

`gemini` will internally use the returned variant_ids to modifiy the sqlite
query to select only those rows.

The `query` function only needs to worry about which variant_ids to return,
not how to integrate with the rest of `gemini`.

bcolz implementation
====================

We have a reference implementation using `bcolz <http://bcolz.blosc.org/>`_.
It can be used from gemini by running::

    gemini bcolz_index $db

This is easily parallelized by specifying a column per process, e.g.:

     gemini bcolz_index $db --cols gt_types

Which can index nearly 9K variants / second for 17 samples in our testing.

It is recommended to only index the columns you'll be using in the 
`--gt-filter`.

Indexing is done only once to create the bcolz index.
After that, add `--use-bcolz` to an existing gemini query command. e.g.::

    gemini query -q "select variant_id, gts.1719PC0016 from variants"  \
        --gt-filter "gt_types.1094PC0012 == HET and gt_types.1719PC0016 == HET and gts.1094PC0012 == 'A/C'" \
        --use-bcolz \
        test/test.query.db 


This query will return identical results with or without using bcolz. It is likely
only faster to use `bcolz` on complex queries that are slow with the default gemini
apparatus.

.. note ::

    indexing the 'gts' column will be much slower (only about 350 variants per second instead of
    up to 25K per second) as it must be stored as an object rather than a fixed-size numeric type.
    It will be a much larger index. So only create an index on 'gts' if necessary.

