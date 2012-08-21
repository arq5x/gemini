#!/usr/bin/env python

import os
import sys
import sqlite3
import numpy as np
import cPickle
import zlib
from collections import defaultdict
from gemini.config import read_gemini_config
import gemini_utils as util

config = read_gemini_config()
path_dirname = config["annotation_dir"]


def get_pathways(args):

    if args.version == None or args.version == '66':
        path_file = os.path.join(path_dirname, 'kegg_pathways_ensembl66') 
    elif args.version == '67':
        path_file = os.path.join(path_dirname, 'kegg_pathways_ensembl67')
    elif args.version == '68':
        path_file = os.path.join(path_dirname, 'kegg_pathways_ensembl68')
    else:
        sys.exit("Unsupported Ensembl gene version.\n")

    agn_paths = defaultdict(list)
    hgnc_paths = defaultdict(list)
    ensembl_paths = defaultdict(list)

    for line in open(path_file, 'r'):
        fields=line.strip().split("\t")
        uniprot = fields[0]
        agn = fields[1]
        hgnc = fields[2]
        ensid = fields[3]
        ens_transcript = fields[4]
        hsa = fields[5]
        path = fields[6] if fields[6] != 'None' else None

        # clean up the pathways such that this:
        # path:hsa00260;Glycine_serine_and_threonine_metabolism
        # becomes this:
        # hsa00260:Glycine_serine_and_threonine_metabolism
        if path is not None and path.startswith("path:"):
            path = path[5:]
            path = path.replace(";", ":")

        # build gene/transcript -> pathway mappings using
        # all three gene naming conventions
        agn_paths[(agn, ens_transcript)].append(path)
        hgnc_paths[(hgnc, ens_transcript)].append(path)
        ensembl_paths[(ensid, ens_transcript)].append(path)

    return agn_paths, hgnc_paths, ensembl_paths

def _get_pathways(gene, transcript, pathways, allow_none=True):
    # get distinct pathways
    pathways = set(pathways)
    if None in pathways:
        # remove "None" if a valid pathway exists.
        if len(pathways) > 1 or allow_none is False:
           pathways.remove(None)
    return pathways

def _report_variant_pathways(c, args, idx_to_sample):
    
    (agn_paths, hgnc_paths, ensembl_paths) = get_pathways(args)
    
    for r in c:
        gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
        gts      = np.array(cPickle.loads(zlib.decompress(r['gts'])))        
        gene     = str(r['gene'])
        trans    = str(r['transcript'])
        
        pathways = []
        if (gene, trans) in agn_paths:
            pathways = _get_pathways(gene, trans, agn_paths[(gene, trans)], 
                            allow_none=False)
        elif (gene, trans) in hgnc_paths:
            pathways = _get_pathways(gene, trans, hgnc_paths[(gene, trans)],
                            allow_none=False)
        elif (gene, trans) in ensembl_paths:
            pathways = _get_pathways(gene, trans, ensembl_paths[(gene, trans)],
                            allow_none=False)
        pathlist = ",".join(pathways)
        for idx, type in enumerate(gt_types):
            if type > 0 and len(pathways) > 0:
                print "\t".join([r['chrom'], str(r['start']), str(r['end']), \
                                 r['ref'], r['alt'], r['most_severe_impact'], \
                                 idx_to_sample[idx], gts[idx], gene, trans, \
                                 pathlist])
    
def get_ind_pathways(c, args):

    idx_to_sample = util.map_indicies_to_samples(c)

    query = "SELECT DISTINCT v.chrom, v.start, v.end, v.ref, v.alt, \
                             v.most_severe_impact, v.gt_types, v.gts, i.gene, \
                             i.transcript \
             FROM variants v, variant_impacts i \
             WHERE v.variant_id = i.variant_id"

    c.execute(query)

    # header
    print '\t'.join(['chrom', 'start', 'end', 'ref', 'alt', \
                     'highest_impact', 'sample', 'genotype', \
                     'gene', 'transcript', 'pathway'])

    _report_variant_pathways(c, args, idx_to_sample)
                            

def get_ind_lof_pathways(c, args):

    idx_to_sample = util.map_indicies_to_samples(c)

    query = "SELECT DISTINCT v.chrom, v.start, v.end, v.ref, v.alt, \
                             v.most_severe_impact, v.gt_types, v.gts, i.gene, \
                             i.transcript \
             FROM variants v, variant_impacts i \
             WHERE v.variant_id = i.variant_id \
             AND i.is_lof='1'"

    c.execute(query)

    # header
    print '\t'.join(['chrom', 'start', 'end', 'ref', 'alt', \
                     'highest_impact', 'sample', 'genotype', \
                     'gene', 'transcript', 'pathway'])

    _report_variant_pathways(c, args, idx_to_sample)



def pathways(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        
        if (not args.lof):
            get_ind_pathways(c, args)
        else:
            get_ind_lof_pathways(c, args)



          