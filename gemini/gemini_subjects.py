#!/usr/bin/env python
import sys
from collections import defaultdict
from compiler import compile

from inheritance import Family
import sqlalchemy as sql
import database
from gemini_constants import *
import GeminiQuery

from functools import wraps

def compile_decorator(f):
    """decorator to automatically compile the eval strings returned from
    the filter methods"""
    @wraps(f)
    def wrapper(*args, **kwargs):
        query_string = f(*args, **kwargs)
        if query_string == "False" or query_string == {"any": "False"}:
            return None
        if not isinstance(query_string, dict):
            return compile(query_string, "<string>", "eval")
        query_dict = query_string
        for k, stmt in query_dict.iteritems():
            query_dict[k] = compile(stmt, "<string>", "eval")
        return query_dict
    return wrapper

def get_phred_query(sample_id, gt_ll, genotype, prefix=" and ", invert=False):
    """Default is to test < where a low value phred-scale is high
    confidence for that genotype
    >>> get_phred_query(2, 22, "het")
    ' and gt_phred_ll_het[1] < 22'

    >>> get_phred_query(2, 22, "het", prefix="")
    'gt_phred_ll_het[1] < 22'

    >>> get_phred_query(2, 22, "het", prefix="", invert=True)
    'gt_phred_ll_het[1] > 22'

    """
    assert genotype in ("het", "homref", "homalt")
    if not gt_ll: return ""

    # they passed in the subject:
    if hasattr(sample_id, "sample_id"):
        sample_id = sample_id.sample_id

    sign = ["<", ">"][int(invert)]
    s = "gt_phred_ll_{genotype}[{sample_id}] {sign} {gt_ll}"\
            .format(sample_id=sample_id-1, genotype=genotype,
                    gt_ll=gt_ll, sign=sign)
    return prefix + s


class Subject(object):

    """
    Describe a single subject in the the samples table.
    """
    def __init__(self, row):
        self._set_fields_from_row(row)

    def __repr__(self):
        return "\t".join(map(str, [self.name, self.paternal_id,
                                   self.maternal_id, self.phenotype]))

    def set_father(self):
        self.father = True

    def set_mother(self):
        self.mother = True

    def _set_fields_from_row(self, row):
        self.__dict__.update(row)
        #for k, v in zip(row.keys(), row):
        #    self.__dict__[k] = v
        self.phenotype = int(self.phenotype) if self._has_phenotype() else None
        self._set_affected_status()

    def _has_phenotype(self):
        if hasattr(self, 'phenotype') and self.phenotype is not None:
            return True

    def _set_affected_status(self):
        # 1 = unaffected
        # 2 = affected
        # 0 or -9 is unknown.
        # http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
        pheno = str(self.phenotype)
        if pheno == "2":
            self.affected = True
        elif pheno == "1":
            self.affected = False
        # distinguish unknown from known to be unaffected.
        else:
            self.affected = None

def get_families(db, selected_families=None):
    """
    Query the samples table to return a list of Family
    objects that each contain all of the Subjects in a Family.
    """
    conn, metadata = database.get_session_metadata(db)

    families_dict = Family.from_cursor(conn)

    # if the user has specified a set of selected families
    # to which the analysis should be restricted, then
    # first sanity check that the family ids they specified are valid.
    if selected_families is not None:
        for family in selected_families.split(','):
            if family not in families_dict:
                sys.exit("ERROR: family \"%s\" is not a valid family_id\n" % family)

    families = []
    for fam in families_dict:
        if selected_families is None or fam in selected_families:
            families.append(families_dict[fam])
    return families

def get_family_dict(args):
    families = defaultdict(list)
    subjects = get_subjects(args)
    for subject in subjects.values():
        families[subject.family_id].append(subject)

    return families

def get_subjects(args, skip_filter=False):
    """
    return a dictionary of subjects, optionally using the
    subjects_query argument to filter them.
    """
    gq = GeminiQuery.GeminiQuery(args.db)

    #query = "SELECT * FROM samples"
    query = ""
    if not skip_filter:
        if hasattr(args, 'sample_filter') and args.sample_filter:
            query += args.sample_filter

    res = gq.metadata.tables["samples"].select().where(sql.text(query)).execute()

    samples_dict = {}
    for row in res:
        subject = Subject(row)
        samples_dict[subject.name] = subject
    return samples_dict

def get_subjects_in_family(args, family):
    subjects = get_subjects(args)
    family_names = [f.name for f in family]
    subject_dict = {}
    for subject in subjects:
        if subject in family_names:
            subject_dict[subject] = subjects[subject]
    return subject_dict

if __name__ == "__main__":
    import doctest
    doctest.testmod()
