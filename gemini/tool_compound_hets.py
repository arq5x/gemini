#!/usr/bin/env python
import sqlite3
import os
import sys
import collections
import re
from copy import copy

import GeminiQuery
import gemini_utils as util
from gemini_constants import *
import gemini_subjects as subjects

class Site(object):
    def __init__(self, row):
        self.row = row
        self.phased = None
        self.gt = None

    def __eq__(self, other):
        return self.row['chrom'] == other.row['chrom'] and \
               self.row['start'] == other.row['start']

    def __repr__(self):
        return ",".join([self.row['chrom'], 
                         str(self.row['start']), 
                         str(self.row['end'])])

    def __hash__(self):
        "hash the site based on chrom+start"
        return sum(ord(c) for c in self.row['chrom']) + int(self.row['start'])


def _add_necessary_columns(args, custom_columns):
    """
    Convenience function to tack on columns that are necessary for
    the functionality of the tool but yet have not been specifically
    requested by the user.
    """
    # we need to add the variant's chrom, start and gene if 
    # not already there.
    if custom_columns.find("gene") < 0:
        custom_columns += ", gene"
    if custom_columns.find("start") < 0:
        custom_columns += ", start"
    if custom_columns.find("alt") < 0:
        custom_columns += ", alt"
        
    return custom_columns

def create_query(args):
    """
    Construct a query to identify candidate compound heterozygotes
    based on the user's columns and filters
    """
    if args.columns is not None:
        custom_columns = _add_necessary_columns(args, str(args.columns))        
        query = "SELECT " + custom_columns + \
                " FROM variants " + \
                " WHERE (is_exonic = 1 or impact_severity != 'LOW') "
    else:
        # report the kitchen sink
        query = "SELECT *" + \
                ", gts, gt_types, gt_phases, gt_depths, \
                gt_ref_depths, gt_alt_depths, gt_quals" + \
                " FROM variants " + \
                " WHERE (is_exonic = 1 or impact_severity != 'LOW') "

    # add any non-genotype column limits to the where clause
    if args.filter:
        query += " AND " + args.filter
    
    # we need to order results by gene so that we can sweep through the results
    query += " ORDER BY gene"
    return query


def find_valid_het_pairs(args, sample_hets, subjects_dict):
    """
    Identify candidate heterozygote pairs.
    """  
    samples_w_hetpair = collections.defaultdict(list)
    for sample in sample_hets:
        for gene in sample_hets[sample]:

            # we only care about combinations, not permutations
            # (e.g. only need site1,site2, not site1,site2 _and site2,site1)
            # thus we can do this in a ~ linear pass instead of a ~ N^2 pass
            for idx, site1 in enumerate(sample_hets[sample][gene]):
                for site2 in sample_hets[sample][gene][idx + 1:]:

                    # expand the genotypes for this sample at each site into 
                    # it's composite alleles.  e.g. A|G -> ['A', 'G']
                    alleles_site1 = []
                    alleles_site2 = []
                    if not args.ignore_phasing:
                        alleles_site1 = site1.gt.split('|')
                        alleles_site2 = site2.gt.split('|')
                    else:
                        # split on phased (|) or unphased (/) genotypes
                        alleles_site1 = re.split('\||/', site1.gt)
                        alleles_site2 = re.split('\||/', site2.gt)
                    
                    # it is only a true compound heterozygote IFF
                    # the alternates are on opposite haplotypes.
                    if not args.ignore_phasing:
                        # return the haplotype on which the alternate allele 
                        # was observed for this sample at each candidate het. 
                        # site. e.g., if ALT=G and alleles_site1=['A', 'G']
                        # then alt_hap_1 = 1.  if ALT=A, then alt_hap_1 = 0
                        if "," in str(site1.row['alt']) or \
                           "," in str(site2.row['alt']):
                            sys.stderr.write("WARNING: Skipping candidate for sample"
                                             " %s b/c variants with mult. alt."
                                             " alleles are not yet supported. The sites are:"
                                             " %s and %s.\n" % (sample, site1, site2))
                            continue

                        alt_hap_1 = alleles_site1.index(site1.row['alt'])
                        alt_hap_2 = alleles_site2.index(site2.row['alt'])

                    # keep as a candidate if 
                    #   1. phasing is considered AND the alt alleles are on
                    #      different haplotypes
                    #   2. the user doesn't care about phasing.
                    #
                    # TODO: Phase based on parental genotypes.
                    if (not args.ignore_phasing and alt_hap_1 != alt_hap_2) \
                        or args.ignore_phasing:
                        samples_w_hetpair[(site1,site2)].append(sample)
    
    # send the remaining candidates off for final filtering
    filter_candidates(args, samples_w_hetpair, subjects_dict)


def filter_candidates(args, samples_w_hetpair, subjects_dict):
    """
    Refine candidate heterozygote pairs based on user's filters.
    """  
    comp_het_id = 1
    for comp_het in samples_w_hetpair:

        # count how many families had this comp_het in *solely the affecteds*
        family_affected = collections.Counter()
        family_unaffected = collections.Counter()
        for s in samples_w_hetpair[comp_het]:
            family =  subjects_dict[s].family_id
            if subjects_dict[s].affected:
                family_affected[family] += 1
            else:
                family_unaffected[family] += 1

        # how many affected and unaffected individuals had the comp_het in total?
        num_affected = sum(family_affected[f] for f in family_affected)
        num_unaffected = sum(family_unaffected[f] for f in family_unaffected)

        # how many distinct families had the compound het at all? 
        num_fam_w_comp_het = sum(family_affected[f] > 0 for f in family_affected) + \
                             sum(family_unaffected[f] > 0 for f in family_unaffected)
        # how many distinct families had the compound het in its affected 
        # family members and NOT in its unaffected?
        num_fam_w_comp_het_strict = 0
        for f in family_affected:
            if family_affected[f] > 0 and \
            (family_unaffected[f] == 0 or family_unaffected[f] is None):
                num_fam_w_comp_het_strict += 1

        if args.only_affected:
            if num_unaffected > 0:
                continue
            if num_fam_w_comp_het_strict < args.min_kindreds:
                continue
        else:
            if num_fam_w_comp_het < args.min_kindreds:
                continue

        for s in samples_w_hetpair[comp_het]:
            family =  subjects_dict[s].family_id
            if args.families is not None and family not in args.families.split(','):
                continue
            else:
                print "\t".join([str(family), s, str(comp_het_id),
                             str(comp_het[0].row)])
                print "\t".join([str(family), s, str(comp_het_id),
                              str(comp_het[1].row)])
        comp_het_id += 1


def get_compound_hets(args):
    """
    Report candidate compound heterozygotes.
    """
    gq = GeminiQuery.GeminiQuery(args.db, include_gt_cols=True)
    idx_to_sample = gq.idx_to_sample
    subjects_dict = subjects.get_subjects(args)
    
    # run the query applying any genotype filters provided by the user.
    gq.run(create_query(args))

    sample_hets = collections.defaultdict(lambda: collections.defaultdict(list))
    curr_gene = None
    prev_gene = None
    # output header
    print "family\tsample\tcomp_het_id\t" + str(gq.header)
    # Collect all of the genic heterozygotes for each sample / gene
    for row in gq:

        gt_types = row['gt_types']
        gt_bases = row['gts']
        gt_phases = row['gt_phases']
        curr_gene = row['gene']
        
        # gene has changed. reset.
        if curr_gene != prev_gene and prev_gene is not None:
            find_valid_het_pairs(args, sample_hets, subjects_dict)
            sample_hets = collections.defaultdict(lambda: collections.defaultdict(list))
        else:        
            site = Site(row)

            # track each sample that is heteroyzgous at this site.
            for idx, gt_type in enumerate(gt_types):
                if gt_type == HET:
                    sample = idx_to_sample[idx]
                    sample_site = copy(site)
                    sample_site.phased = gt_phases[idx]

                    if not sample_site.phased and not args.ignore_phasing:
                        continue

                    sample_site.gt = gt_bases[idx]
                    # add the site to the list of candidates for this sample/gene
                    sample_hets[sample][site.row['gene']].append(sample_site)

        prev_gene = curr_gene

    # process the last gene seen
    find_valid_het_pairs(args, sample_hets, subjects_dict)    


def run(parser, args):
    if os.path.exists(args.db):
        get_compound_hets(args)
        
