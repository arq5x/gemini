import sys
import os
from compiler import compile
import operator
import itertools as it
from argparse import ArgumentParser
from GeminiQuery import GeminiQuery


def add_args(a=None):
    if a is None:
        a = ArgumentParser()
    a.add_argument("--min-filters", type=int, default=1)
    a.add_argument("--gt-filter", required=True, default=[], action='append')
    a.add_argument("--filter")
    a.add_argument("--columns", default="chrom,start,end,gene,impact,impact_severity,max_aaf_all")
    a.add_argument("db")
    return a


def add_cols(cols, gt_filter):
    assert isinstance(cols, list)
    all_cols = ["gts", "gt_types", "gt_phases", "gt_depths",
                "gt_ref_depths", "gt_alt_depths", "gt_quals",
                "gt_phred_ll_homref", "gt_phred_ll_het",
                "gt_phred_ll_homalt"]
    return [x for x in all_cols if x in gt_filter and not x in cols]


def gen_results(rows, gt_filters, min_filters, min_variants, columns,
        user_dict=None):
    # we track the index of the passed filter in passed_filters.
    gene_passed_filters = {}
    if user_dict is None:
        user_dict = {}
    subset = []
    for row in rows:
        cols = {c: row[c] for c in columns}
        cols.update(user_dict)
        # have to test all filters since 1 row can meet multiple fitlers.
        row_passed_filters = []
        for i, gt_filter in enumerate(gt_filters, start=1):
            if eval(gt_filter, cols):
                # track that this filter passed.
                gene_passed_filters[i] = True
                row_passed_filters.append(i)
        if row_passed_filters:
            row.print_fields['variant_filters'] = ",".join(map(str, row_passed_filters))
            subset.append(row)
    if len(gene_passed_filters) < min_filters or len(subset) < min_variants:
        raise StopIteration

    # e.g. 1,2 indicating which filters passed
    passed_filters = ",".join(str(x) for x in sorted(gene_passed_filters))
    for row in subset:
        row.print_fields['n_gene_variants'] = len(subset)
        row.print_fields['gene_filters'] = passed_filters
        yield row


def genewise(db, gt_filters, filter=None, columns=None, min_filters=None,
             min_variants=1,
             grouper="gene"):
    assert os.path.exists(db)

    orig_columns = [x.strip() for x in (columns or "").split(",")]
    added_cols = add_cols(orig_columns, "||".join(gt_filters))
    if grouper not in orig_columns:
        added_cols.append(grouper)
    columns = orig_columns + added_cols
    assert not any(';' in c for c in columns)

    # NOTE: we could make the WHERE part customizable.
    query = "SELECT {columns} FROM variants WHERE (is_exonic = 1 AND impact_severity != 'LOW')"
    if filter:
        query += " AND  " + filter
    query += " ORDER BY CHROM, gene"

    gq = GeminiQuery(db, include_gt_cols=True)

    # use the GeminiQuery machinery to transform to something that's eval'able
    # then compile once for speed.
    cleaned_filters = []
    for gt_filter in gt_filters:
        assert gq._is_gt_filter_safe(gt_filter)
        gt_filter = gq._correct_genotype_filter(gt_filter)
        cleaned_filters.append(compile(gt_filter, gt_filter, 'eval'))

    gq.run(query.format(columns=", ".join(columns)))

    if isinstance(grouper, basestring):
        grouper = operator.itemgetter(grouper)

    user_dict = dict(sample_info=gq.sample_info)
    header_printed = False
    for groupkey, grp in it.groupby(gq, grouper):
        grp = list(grp)
        for x in gen_results(list(grp), cleaned_filters, min_filters or 0,
                             min_variants, columns, user_dict=user_dict):
            for c in added_cols:
                if c != 'gene':
                    del x.print_fields[c]
            if not header_printed:
                print "\t".join(x.print_fields.keys())
                header_printed = True
            print x


def run(args):
    if args.min_filters > len(args.gt_filter):
        sys.stderr.write("ERROR gene-wise: specified --min-filter > the number of --gt-filters\n")
        sys.exit(2)
    genewise(args.db, args.gt_filter, args.filter, args.columns,
             args.min_filters)
