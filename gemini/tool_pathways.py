#!/usr/bin/env python

import os
import sys
import sqlite3
from collections import defaultdict
from gemini.config import read_gemini_config

config = read_gemini_config()
path_dirname = config["annotation_dir"]

def get_pathways(args):     
        
    if args.version == None or args.version == 'v66':
        path_file = os.path.join(path_dirname, 'kegg_pathways_ensembl66') 
    elif args.version == 'v67':
        path_file = os.path.join(path_dirname, 'kegg_pathways_ensembl67')
    elif args.version == 'v68':
        path_file = os.path.join(path_dirname, 'kegg_pathways_ensembl68')
    else:
        sys.stderr.write("Unsupported version\n")
    
    if args.version in ('v66', 'v67', 'v68'):       
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
            path = fields[6]
            # build gene/transcript -> pathway mappings using
            # all three gene naming conventions
            agn_paths[(agn, ens_transcript)].append(path)
            hgnc_paths[(hgnc, ens_transcript)].append(path)
            ensembl_paths[(ensid, ens_transcript)].append(path)
            
        return agn_paths, hgnc_paths, ensembl_paths

def _print_pathways(gene, transcript, pathways):
    for pathway in pathways:
        if pathway is not None:
            print "\t".join([gene, transcript, pathway])
        else:
            print "\t".join([gene, transcript, 'None'])

def get_all_pathways(c, args):
    
    (agn_paths, hgnc_paths, ensembl_paths) = get_pathways(args)

    query = "SELECT DISTINCT gene, transcript from variant_impacts"
    c.execute(query)

    for row in c:
        gene = str(row['gene'])
        trans = str(row['transcript'])

        if (gene, trans) in agn_paths:
            _print_pathways(gene, trans, agn_paths[(gene, trans)])
        elif (gene, trans) in hgnc_paths:
            _print_pathways(gene, trans, hgnc_paths[(gene, trans)])
        elif (gene, trans) in ensembl_paths:
            _print_pathways(gene, trans, ensembl_paths[(gene, trans)])
        else:
            _print_pathways(gene, trans, [None])


def get_lof_pathways(c, args):
    
    (agn_paths, hgnc_paths, ensembl_paths) = get_pathways(args)

    query = "SELECT DISTINCT gene, transcript \
             FROM variant_impacts where is_lof='1'"
    c.execute(query)
    
    for row in c:
        gene = str(row['gene'])
        trans = str(row['transcript'])
    
        if (gene, trans) in agn_paths:
            _print_pathways(gene, trans, agn_paths[(gene, trans)])
        elif (gene, trans) in hgnc_paths:
            _print_pathways(gene, trans, hgnc_paths[(gene, trans)])
        elif (gene, trans) in ensembl_paths:
            _print_pathways(gene, trans, ensembl_paths[(gene, trans)])   


def allquery(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        get_all_pathways(c, args)

def lofquery(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        get_lof_pathways(c, args)

          