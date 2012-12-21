#!/usr/bin/env python
import sqlite3
import os
import sys
import numpy as np
import cPickle
import zlib
import collections
from copy import copy

import gemini_utils as util
from gemini_constants import *

class Subject(object):
    def __init__(self, row):
        self.sample_id = row['sample_id']
        self.name = row['name']
        self.family_id = row['family_id']
        self.paternal_id = row['paternal_id']
        self.maternal_id = row['maternal_id']
        self.sex = row['sex']
        self.phenotype = int(row['phenotype'])
        self.ethnicity = row['ethnicity']
        
        # 1 = unaffected
        # 2 = affected
        # 0 or -9 is unknown.
        # http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
        if self.phenotype == 2:
            self.affected = True
        else:
            self.affected = False
        
    def __repr__(self):
        return "\t".join([self.name, str(self.phenotype)])
        
    def set_father(self):
        self.father = True
        
    def set_mother(self):
        self.mother = True
        
class Family(object):
    def __init__(self, subjects):
        self.subjects = subjects
        self.father = None
        self.mother = None
        self.family_id = self.subjects[0].family_id
        self.children = []
        
    def find_parents(self):
        """
        Screen for children with parental ids so that
        we can identify the parents in this family.
        
        NOTE: assumes at most a 2 generation family.
        """
        for subject in self.subjects:
            # if mom and dad are found, we know this is the child
            if subject.maternal_id is not None and \
               subject.maternal_id != -9 and \
               subject.paternal_id is not None and \
               subject.paternal_id != -9:
                self.father_name = str(subject.paternal_id)
                self.mother_name = str(subject.maternal_id)
                self.children.append(subject)
        
        # now track the actual sampleIds for the parents
        for subject in self.subjects:
            if subject.name == self.father_name:
                self.father = subject
            elif subject.name == self.mother_name:
                self.mother = subject

    def get_auto_recessive_filter(self):
        """
        Generate an eval() filter to apply for this family.
        For example:
        # mom                dad                 affected kid
        '(gt_types[57] == 1 and gt_types[58] == 1 and gt_types[11] == 3)'
        """
        
        # identify which samples are the parents in the family.
        self.find_parents()
        
        mask = "("
        mask += 'gt_types[' + str(self.father.sample_id - 1) + "] == " + \
                 str(GT_HET)
        mask += " and "
        mask += 'gt_types[' + str(self.mother.sample_id - 1) + "] == " + \
                 str(GT_HET)
        mask += " and "
        for i, child in enumerate(self.children):
            if child.affected:
                mask += 'gt_types[' + str(child.sample_id - 1) + "] == " + \
                         str(GT_HOM_ALT)
            else:
                mask += 'gt_types[' + str(child.sample_id - 1) + "] != " + \
                         str(GT_HOM_ALT)
            
            if i < (len(self.children) - 1):
                mask += " and "
        
        mask += ")"
        return mask
        
    def get_subject_columns(self):
        columns = []
        columns.append('gts[' + str(self.father.sample_id - 1) + ']')
        columns.append('gts[' + str(self.mother.sample_id - 1) + ']')
        for child in self.children:
            columns.append('gts[' + str(child.sample_id - 1) + ']')
        
        return columns
        
    def get_subject_labels(self):
        subjects = []
        subjects.append(self.father.name + "(father)")
        subjects.append(self.mother.name + "(mother)")
        for child in self.children:
            if child.affected is True:
                subjects.append(child.name + "(child; affected)")
            else:
                subjects.append(child.name + "(child; unaffected)")
        
        return subjects


def get_families(c):
    query = "SELECT * FROM samples \
             WHERE family_id is not NULL \
             ORDER BY family_id"
    c.execute(query)
    
    families_dict = {}
    for row in c:
        subject = Subject(row)
        family_id = subject.family_id
        if family_id in families_dict:
            families_dict[family_id].append(subject)
        else:
            families_dict[family_id] = []
            families_dict[family_id].append(subject)
    
    families = []
    for fam in families_dict:
        family = Family(families_dict[fam])
        families.append(family)
    return families
    

def get_auto_recessive_candidates(c, families):
    """
    Report candidate compound heterozygous mutations.
    """
    
    for family in families:
        
        query = "SELECT chrom, start, end, ref, alt, gene, \
                        impact, impact_severity, gt_types, gts \
                 FROM variants \
                 WHERE impact_severity != 'LOW'"

        c.execute(query)
        all_query_cols = [str(tuple[0]) for tuple in c.description \
                                            if not tuple[0].startswith("gt")]
                                  
        family_genotype_mask  = family.get_auto_recessive_filter()
        family_sample_columns = family.get_subject_columns()
        family_sample_labels = family.get_subject_labels()
        
        # print a header
        print "=========================="
        print "FAMILY:", family.family_id
        print "=========================="
        print '\t'.join(col for col in all_query_cols),
        print '\t'.join(col for col in family_sample_labels)
        
        # report the resulting auto_rec variants for this familiy
        for row in c:
            gt_types  = \
                np.array(cPickle.loads(zlib.decompress(row['gt_types'])))
            gts  = \
                np.array(cPickle.loads(zlib.decompress(row['gts'])))
                
            # does the variant meet the inheritance model for this family?
            if not eval(family_genotype_mask):
                continue
        
            for col in all_query_cols:
                if col == 'gt_types' or col == 'gts':
                    continue
                print str(row[col]) + '\t',
            for col in family_sample_columns:
                print str(eval(col)) + '\t',
            print
                
                


def run(parser, args):

    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        families = get_families(c)
        get_auto_recessive_candidates(c, families)



