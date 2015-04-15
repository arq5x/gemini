#!/usr/bin/env python
import sqlite3
import os
import sys
import numpy as np
from collections import defaultdict
import itertools as it
from compiler import compile

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
        self.__dict__.update(dict(zip(row.keys(), row)))
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


class Family(object):

    """
    Describe the relationships among multiple subjects in a family.
    """
    def __init__(self, subjects):
        self.subjects = subjects
        self.father = None
        self.mother = None
        self.family_id = self.subjects[0].family_id
        self.children = []
        self.affected = []
        self.unaffected = []
        self.affected_children = []
        self.unaffected_children = []

        self.is_constructed = False

        self.build_family()


    def has_an_affected(self):
        """
        Return True if the Family has at least one affected individual.
        Otherwise return False.
        """
        for subject in self.subjects:
            if subject.affected:
                return True
        return False

    def has_an_affected_child(self):
        """
        Return True if the Family has at least one affected child.
        Otherwise return False.
        """
        if not self.is_constructed:
            self.build_family()

        for child in self.children:
            if child.affected:
                return True
        return False

    def build_family(self):
        """
        Screen for children with parental ids so that
        we can identify the parents in this family.

        NOTE: assumes at most a 2 generation family.
        """

        # build only if the family has not already been built.
        if self.is_constructed is False:
            self.father_name = None
            self.mother_name = None
            for subject in self.subjects:
                # if mom and dad are found, we know this is the child
                if subject.maternal_id != "0" and subject.paternal_id != "0":
                    self.father_name = str(subject.paternal_id)
                    self.mother_name = str(subject.maternal_id)
                    self.children.append(subject)

            # now track the actual sampleIds for the parents
            for subject in self.subjects:
                if self.father_name is not None and \
                   subject.name == self.father_name:
                    self.father = subject
                elif self.mother_name is not None and \
                   subject.name == self.mother_name:
                    self.mother = subject

            # prevent reconstructing family every time function is called.
            self.is_constructed = True

        if self.father is not None and self.mother is not None:
            return True
        else:
            return False

    @compile_decorator
    def get_auto_recessive_filter(self, gt_ll=False):
        """
        Generate an autosomal recessive eval() filter to apply for this family.
        For example:

        '(gt_types[57] == HET and \  # mom
          gt_types[58] == HET and \  # dad
          gt_types[11] == HOM_ALT)'  # affected child
        """

        parents_found = self.build_family()
        affected_found = self.has_an_affected()
        # identify which samples are the parents in the family.
        # Fail if both parents are not found
        if not parents_found and not affected_found:
            sys.stderr.write("WARNING: Unable to identify at least one "
                             "affected individual for family (%s). "
                             "Consequently, GEMINI will not screen for "
                             "variants in this family.\n"
                 % self.family_id)
            return "False"

        elif not parents_found and affected_found:
            sys.stderr.write("WARNING: Unable to identify parents for family (%s). "
                             "Consequently, GEMINI will solely place genotype "
                             "requirements on subjects based on their phenotype.\n"
                             % self.family_id)

            mask = "("
            for i, subject in enumerate(self.subjects):
                if subject.affected:
                    mask += '(gt_types[' + str(subject.sample_id - 1) + "] == " + \
                        str(HOM_ALT)
                    mask += get_phred_query(subject, gt_ll, "homalt") + ")"
                else:
                    mask += '(gt_types[' + str(subject.sample_id - 1) + "] != " + \
                        str(HOM_ALT)
                    mask += get_phred_query(subject, gt_ll, "homalt", invert=True) + ")"

                if i < (len(self.subjects) - 1):
                    mask += " and "
            mask += ")"

            return mask

        elif parents_found:
            # if either parent is affected, this family cannot satisfy
            # a recessive model, as the parents should be carriers.
            if self.father.affected is True or self.mother.affected is True:
                return "False"

            mask = "("
            mask += 'gt_types[' + str(self.father.sample_id - 1) + "] == " + \
                str(HET)

            mask += get_phred_query(self.father.sample_id, gt_ll, "het")

            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + \
                str(HET)

            mask += get_phred_query(self.mother.sample_id, gt_ll, "het")

            if self.has_an_affected_child():
                for i, child in enumerate(self.children):
                    if child.affected is True:
                        mask += " and "
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + \
                            str(HOM_ALT)

                        mask += get_phred_query(child.sample_id, gt_ll, "homalt")

                    # only allow an unaffected if there are other affected children
                    elif child.affected is False and self.has_an_affected_child():
                        mask += " and "
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] != " + \
                            str(HOM_ALT)

                        # TODO: # Aaron, does this seem right? We want them not
                        # HOMALT, so we want a high GT_LL for the HOM ALT
                        # allele.
                        mask += get_phred_query(child.sample_id, gt_ll, "homalt", invert=True)
                    elif child.affected is None:
                        # assume just testing for inheritance patterns
                        mask += " and "
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + \
                            str(HOM_ALT)
                        mask += get_phred_query(child.sample_id, gt_ll, "homalt")
            else:
                for i, child in enumerate(self.children):
                    mask += " and "
                    mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + \
                        str(HOM_ALT)
                    mask += get_phred_query(child.sample_id, gt_ll, "homalt")

            mask += ")"

            return mask

    @compile_decorator
    def get_auto_dominant_filter(self, gt_ll=False):
        """
        Generate an autosomal dominant eval() filter to apply for this family.
        For example:
        '(
          ((bool(gt_types[57] == HET)         # mom
            != \
            bool(gt_types[58] == HET)) and \  # dad
            gt_types[11] == HET               # affected child
        )'

        NOTE: the bool(dad) != bool(mom) is an XOR requiring that one and
        only one of the parents is heterozygous
        """

        parents_found = self.build_family()
        affected_found = self.has_an_affected()

        # identify which samples are the parents in the family.
        # Fail if both parents are not found
        if not parents_found and not affected_found:
            sys.stderr.write("WARNING: Unable to identify at least one "
                             "affected individual for family (%s). "
                             "Consequently, GEMINI will not screen for "
                             "variants in this family.\n"
                 % self.family_id)
            return "False"

        elif not parents_found and affected_found:
            sys.stderr.write("WARNING: Unable to identify parents for family (%s). "
                             "Consequently, GEMINI will solely place genotype "
                             "requirements on subjects based on their phenotype.\n"
                             % self.family_id)

            mask = "("
            for i, subject in enumerate(self.subjects):
                if subject.affected:
                    mask += 'gt_types[' + str(subject.sample_id - 1) + "] == " + str(HET)
                    mask += get_phred_query(subject.sample_id, gt_ll, "het")
                else:
                    mask += 'gt_types[' + str(subject.sample_id - 1) + "] == " + str(HOM_REF)
                    mask += get_phred_query(subject.sample_id, gt_ll, "homref")

                if i < (len(self.subjects) - 1):
                    mask += " and "
            mask += ")"

            return mask

        elif parents_found:
            mask = ""
            if self.father.affected is True and self.mother.affected is True:
                # doesn't meet an auto. dominant model if both parents are affected
                # [*]---(*)
                #     |
                #    (*)
                return "False"
            elif ((self.father.affected is False and self.mother.affected is False)
                 or
                 (self.father.affected is None and self.mother.affected is None)):
                # if neither parents are affected, or the affection status is
                # unknown for both, we can just screen for variants where one and
                # only one of the parents are hets and and the child is also a het
                # []---()
                #    |
                #   (*)
                mask = "((bool("
                mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HET)
                # TODO: ask Aaron. adding these likelihoods means we are
                # checking that only one of them is *confidently* a het.
                # where as without the likelihoods, it could be more stringent..
                mask += get_phred_query(self.father.sample_id, gt_ll, "het") + ")"
                mask += ") != "
                mask += 'bool((gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HET)
                mask += get_phred_query(self.mother.sample_id, gt_ll, "het") + ")"
                mask += ")) and "


                for i, child in enumerate(self.children):
                    if child.affected:
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)
                        mask += get_phred_query(child.sample_id, gt_ll, "het")

                    else:
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + \
                            str(HOM_REF)
                        mask += get_phred_query(child.sample_id, gt_ll, "homref")

                    if i < (len(self.children) - 1):
                        mask += " and "
                mask += ")"
                return mask
            elif (self.father.affected is True and
                  self.mother.affected is not True):
                # if only Dad is known to be affected, we must enforce
                # that only the affected child and Dad have the
                # same heterozygous genotype.
                # [*]---()
                #     |
                #    (*)
                mask = "(("
                mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HET)
                mask += get_phred_query(self.father.sample_id, gt_ll, "het") + ")"
                mask += " and "
                mask += '(gt_types[' + str(self.mother.sample_id - 1) + "] != " + str(HET)
                mask += get_phred_query(self.mother.sample_id, gt_ll, "het", invert=True) + ")"

                mask += ") and "
                for i, child in enumerate(self.children):
                    if child.affected:
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)
                        mask += get_phred_query(child.sample_id, gt_ll, "het")
                    else:
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_REF)
                        mask += get_phred_query(child.sample_id, gt_ll, "homref")
                    if i < (len(self.children) - 1):
                        mask += " and "
                mask += ")"
                return mask
            elif (self.father.affected is not True
                  and self.mother.affected is True):
                # if only Mom is known to be affected, we must enforce
                # that only the affected child and Mom have the
                # same heterozygous genotype.
                # []---(*)
                #    |
                #   (*)
                mask = "(("
                mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HET)
                mask += get_phred_query(self.mother.sample_id, gt_ll, "het")
                mask += " and "
                mask += 'gt_types[' + str(self.father.sample_id - 1) + "] != " + str(HET)
                mask += get_phred_query(self.father.sample_id, gt_ll, "het", invert=True)
                mask += ") and "
                for i, child in enumerate(self.children):
                    if child.affected:
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)
                        mask += get_phred_query(child.sample_id, gt_ll, "het")
                    else:
                        mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + \
                              str(HOM_REF)
                        mask += get_phred_query(child.sample_id, gt_ll, "homref")
                    if i < (len(self.children) - 1):
                        mask += " and "
                mask += ")"
                return mask

    @compile_decorator
    def get_de_novo_filter(self, only_affected=False, gt_ll=False):
        """
        Generate aa de novo mutation eval() filter to apply for this family.
        For example:

        '(gt_types[57] == HOM_REF and \  # mom
          gt_types[58] == HOM_REF and \  # dad
          gt_types[11] == HET)'          # affected child

          # [G/G]---(G/G)
          #       |
          #     (A/G)
        """

        # identify which samples are the parents in the family.
        # Fail if both parents are not found
        if not self.build_family():
            sys.stderr.write("WARNING: Unable to find parents for family (%s). "
                 "GEMINI is currently only able to identify candidates "
                 "from two generational families.\n"
                 % self.family_id)
            return "False"

        mask = "("

        mask += "("
        mask += 'gt_types[' + str(self.father.sample_id - 1) + "] == " + \
            str(HOM_REF)
        mask += " and "
        mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + \
            str(HOM_REF)

        mask += get_phred_query(self.father.sample_id, gt_ll, "homref")
        mask += get_phred_query(self.mother.sample_id, gt_ll, "homref")

        mask += ")"


        mask += " or "

        mask += "("
        mask += 'gt_types[' + str(self.father.sample_id - 1) + "] == " + \
            str(HOM_ALT)
        mask += " and "
        mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + \
            str(HOM_ALT)

        mask += get_phred_query(self.father.sample_id, gt_ll, "homalt")
        mask += get_phred_query(self.mother.sample_id, gt_ll, "homalt")

        mask += ")"

        mask += ")"

        mask += " and ("

        if len(self.children) == 1:
            if only_affected == False or \
            (only_affected == True and self.children[0].affected == True):
                mask += '( gt_types[' + str(self.children[0].sample_id - 1) + "] == " + str(HET)
                mask += get_phred_query(self.children[0].sample_id, gt_ll, "het")
                mask += ")"
        else:
            if only_affected == False:
                for i, child in enumerate(self.children):
                    mask += '( gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)
                    mask += get_phred_query(child, gt_ll, "het")
                    mask += ")"

                    if i < (len(self.children) - 1):
                        mask += " or "
            else:
                # one or more of the affecteds must be HET
                num_affected = sum(child.affected for child in self.children)
                affected = 0
                for child in self.children:
                    if child.affected == True:
                        affected += 1
                        mask += '(gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)
                        mask += get_phred_query(child, gt_ll, "het") + ")"

                        if affected < num_affected:
                            mask += " or "

                mask += ") and ("

                # AND, none of the unaffecteds can be HET
                num_unaffected = sum(not child.affected for child in self.children)
                unaffected = 0
                for i, child in enumerate(self.children):
                    if child.affected is False:
                        unaffected += 1
                        mask += '(gt_types[' + str(child.sample_id - 1) + "] != " + str(HET)
                        mask += get_phred_query(child, gt_ll, "het", invert=True) + ")"
                        if unaffected < num_unaffected:
                            mask += " and "
        mask += ")"
        return mask


    @compile_decorator
    def get_mendelian_violation_filter(self, gt_ll=False):
        """
        Generate Mendelian violation eval() filter to apply for this family.
        For example:

        '(gt_types[57] == HOM_REF and \  # mom
          gt_types[58] == HOM_REF and \  # dad
          gt_types[11] == HET)'          # affected child

          # [G/G]---(G/G)
          #       |
          #     (A/G)
        """

        # identify which samples are the parents in the family.
        # Fail if both parents are not found
        if not self.build_family():
            sys.stderr.write("WARNING: Unable to find parents for family (%s). "
                 "GEMINI is currently only able to identify candidates "
                 "from two generational families.\n"
                 % self.family_id)
            return {"any": "False"}

        # outer start paren
        #masks is keys by the type of violation and values are the filters.
        masks = {}

        for i, child in enumerate(self.children):

            ##################################################
            # Plausible de novos
            ##################################################
            # DAD = HOM_REF; MOM = HOM_REF; KID = HET (De novo)
            mask = '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)

            mask += get_phred_query(self.father, gt_ll, "homref")
            mask += get_phred_query(self.mother, gt_ll, "homref")
            mask += get_phred_query(child, gt_ll, "het")


            mask += ")"
            mask += " or "

            # DAD = HOM_ALT; MOM = HOM_ALT; KID = HET (De novo)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HET)

            mask += get_phred_query(self.father, gt_ll, "homalt")
            mask += get_phred_query(self.mother, gt_ll, "homalt")
            mask += get_phred_query(child, gt_ll, "het")

            mask += ")"
            masks['plausible de novo'] = mask

            ##################################################
            # Implausible de novos
            ##################################################

            # DAD = HOM_REF; MOM = HOM_REF; KID = HOM_ALT (Implausible de novo)
            mask = '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_ALT)

            mask += get_phred_query(self.father, gt_ll, "homref")
            mask += get_phred_query(self.mother, gt_ll, "homref")
            mask += get_phred_query(child, gt_ll, "homalt")

            mask += ")"
            mask += " or "

            # DAD = HOM_ALT; MOM = HOM_ALT; KID = HOM_REF (Implausible de novo)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_REF)

            mask += get_phred_query(self.father, gt_ll, "homalt")
            mask += get_phred_query(self.mother, gt_ll, "homalt")
            mask += get_phred_query(child, gt_ll, "homref")

            mask += ")"

            masks['implausible de novo'] = mask

            ##################################################
            # Uniparental disomies
            ##################################################

            # DAD = HOM_REF; MOM = HOM_ALT; KID = HOM_REF (Uniparental disomy)
            mask = '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_REF)

            mask += get_phred_query(self.father, gt_ll, "homref")
            mask += get_phred_query(self.mother, gt_ll, "homalt")
            mask += get_phred_query(child, gt_ll, "homref")

            mask += ")"
            mask += " or "

            # DAD = HOM_REF; MOM = HOM_ALT; KID = HOM_ALT (Uniparental disomy)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_ALT)

            mask += get_phred_query(self.father, gt_ll, "homref")
            mask += get_phred_query(self.mother, gt_ll, "homalt")
            mask += get_phred_query(child, gt_ll, "homalt")

            mask += ")"
            mask += " or "

            # DAD = HOM_ALT; MOM = HOM_REF; KID = HOM_REF (Uniparental disomy)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_REF)

            mask += get_phred_query(self.father, gt_ll, "homalt")
            mask += get_phred_query(self.mother, gt_ll, "homref")
            mask += get_phred_query(child, gt_ll, "homref")

            mask += ")"
            mask += " or "

            # DAD = HOM_ALT; MOM = HOM_REF; KID = HOM_ALT (Uniparental disomy)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_ALT)

            mask += get_phred_query(self.father, gt_ll, "homalt")
            mask += get_phred_query(self.mother, gt_ll, "homref")
            mask += get_phred_query(child, gt_ll, "homalt")

            mask += ")"
            masks['uniparental disomy'] = mask
            ##################################################
            # Losses of heterozygosity
            ##################################################

            # DAD = HOM_REF; MOM = HET; KID = HOM_ALT (Loss of heterozygosity)
            mask = '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HET)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_ALT)

            mask += get_phred_query(self.father, gt_ll, "homref")
            mask += get_phred_query(self.mother, gt_ll, "het")
            mask += get_phred_query(child, gt_ll, "homalt")

            mask += ")"
            mask += " or "

            # DAD = HOM_ALT; MOM = HET; KID = HOM_REF (Loss of heterozygosity)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HET)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_REF)

            mask += get_phred_query(self.father, gt_ll, "homalt")
            mask += get_phred_query(self.mother, gt_ll, "het")
            mask += get_phred_query(child, gt_ll, "homref")

            mask += ")"
            mask += " or "

            # DAD = HET; MOM = HOM_REF; KID = HOM_ALT (Loss of heterozygosity)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HET)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_REF)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_ALT)

            mask += get_phred_query(self.father, gt_ll, "het")
            mask += get_phred_query(self.mother, gt_ll, "homref")
            mask += get_phred_query(child, gt_ll, "homalt")

            mask += ")"
            mask += " or "

            # DAD = HET; MOM = HOM_ALT; KID = HOM_REF (Loss of heterozygosity)
            mask += '(gt_types[' + str(self.father.sample_id - 1) + "] == " + str(HET)
            mask += " and "
            mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + str(HOM_ALT)
            mask += " and "
            mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + str(HOM_REF)

            mask += get_phred_query(self.father, gt_ll, "het")
            mask += get_phred_query(self.mother, gt_ll, "homalt")
            mask += get_phred_query(child, gt_ll, "homref")

            mask += ")"
            masks['loss of heterozygosity'] = mask
        # outer end paren
        return masks


    def get_gt_field(self, field):
        columns = []
        prefix = field + "["
        if not self.build_family():
            for subject in self.subjects:
                columns.append(prefix + str(subject.sample_id - 1) + ']')
        else:
            columns.append(prefix + str(self.father.sample_id - 1) + ']')
            columns.append(prefix + str(self.mother.sample_id - 1) + ']')
            for child in self.children:
                columns.append(prefix + str(child.sample_id - 1) + ']')
        return columns

    def get_genotype_columns(self):
        return self.get_gt_field("gts")

    def get_genotype_depths(self):
        return self.get_gt_field("gt_depths")

    def get_genotype_lls(self):
        return dict(
            homref=self.get_gt_field("gt_phred_ll_homref"),
            het=self.get_gt_field("gt_phred_ll_het"),
            homalt=self.get_gt_field("gt_phred_ll_homalt"))

    def get_genotype_labels(self):
        """
        Return header genotype labels for the parents and the children.
        """
        labels = []

        # these are just anonymous affected and unaffected i
        # individuals in the same family
        if not self.build_family():
            for subject in self.subjects:
                if subject.affected is True:
                    labels.append(subject.name + "(affected)")
                elif subject.affected is False:
                    labels.append(subject.name + "(unaffected)")
                elif subject.affected is None:
                    labels.append(subject.name + "(unknown)")
        else:

            if self.father.affected is True:
                labels.append(self.father.name + "(dad;affected)")
            elif self.father.affected is False:
                labels.append(self.father.name + "(dad;unaffected)")
            elif self.father.affected is None:
                labels.append(self.father.name + "(dad;unknown)")

            if self.mother.affected is True:
                labels.append(self.mother.name + "(mom;affected)")
            elif self.mother.affected is False:
                labels.append(self.mother.name + "(mom;unaffected)")
            elif self.mother.affected is None:
                labels.append(self.mother.name + "(mom;unknown)")

            # handle the childrem

            for child in self.children:
                if child.affected is True:
                    labels.append(child.name + "(child;affected)")
                elif child.affected is False:
                    labels.append(child.name + "(child;unaffected)")
                elif child.affected is None:
                    labels.append(child.name + "(child;unknown)")

        return labels

    def get_subject_depth_labels(self):
        """
        Return header depth labels for the parents and the children.
        """
        subjects = []
        subjects.append(self.father.name + "(depth)")
        subjects.append(self.mother.name + "(depth)")
        for child in self.children:
            subjects.append(child.name + "(depth)")

        return subjects


def get_families(db, selected_families=None):
    """
    Query the samples table to return a list of Family
    objects that each contain all of the Subjects in a Family.
    """
    conn = sqlite3.connect(db)
    conn.isolation_level = None
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    query = "SELECT * FROM samples \
             WHERE family_id is not NULL \
             ORDER BY family_id"
    c.execute(query)

    # create a mapping of family_id to the list of
    # individuals that are members of the family.
    families_dict = {}
    for row in c:
        subject = Subject(row)
        family_id = subject.family_id
        if family_id in families_dict:
            families_dict[family_id].append(subject)
        else:
            families_dict[family_id] = [subject]

    # if the user has specified a set of selected families
    # to which the analysis should be restricted, then
    # first sanity check that the family ids they specified are valid.
    if selected_families is not None:
        for family in selected_families.split(','):
            if family not in families_dict:
                sys.exit("ERROR: family \"%s\" is not a valid family_id\n" % family)

    families = []
    for fam in families_dict:
        if selected_families is None:
            family = Family(families_dict[fam])
            families.append(family)
        elif fam in selected_families:
            family = Family(families_dict[fam])
            families.append(family)
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
    query = "SELECT * FROM samples"
    if not skip_filter:
        if hasattr(args, 'sample_filter') and args.sample_filter:
            query += " WHERE " + args.sample_filter
    gq.c.execute(query)
    samples_dict = {}
    for row in gq.c:
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
