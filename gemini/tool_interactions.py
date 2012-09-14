#!/usr/bin/env python

###########################################################################################################################################
#1. user defined root(may or may not be a mutated gene for the sample); get all interacting partners from the variant list of the sample
#2. For a root fixed to a lof gene in each sample, get the interacting partners
#3. a network graph of mutated genes for each sample, a subnetwork for any gene in that list (not executed here)
###########################################################################################################################################

import os
import sys
import sqlite3
import numpy as np
import zlib
import cPickle
from gemini.config import read_gemini_config
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
from pygraph.algorithms.searching import breadth_first_search
from pygraph.algorithms.minmax import shortest_path
from pygraph.classes.exceptions import AdditionError
from pygraph.algorithms.filters.radius import radius
from pygraph.classes.digraph import digraph
import gemini_utils as util
from gemini_constants import *
from collections import defaultdict

config = read_gemini_config()
path_dirname = config["annotation_dir"]
out = open("file.dot", 'w')
sam = defaultdict(list)
lof = defaultdict(list)

def get_variant_genes(c, args, idx_to_sample):
    for r in c:
        gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
        gts      = np.array(cPickle.loads(zlib.decompress(r['gts'])))
        var_id = str(r['variant_id'])   
        gene     = str(r['gene'])
        impact = str(r['impact'])
        biotype = str(r['biotype'])
        
        for idx, gt_type in enumerate(gt_types):
            if (gt_type == GT_HET or gt_type == GT_HOM_ALT):
                if gene != "None":
                    (key, value) = (idx_to_sample[idx], (gene,var_id,impact,biotype))
                    sam[idx_to_sample[idx]].append(value)
    return sam
    
def get_lof_genes(c, args, idx_to_sample):
    for r in c:
        gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
        gts      = np.array(cPickle.loads(zlib.decompress(r['gts'])))  
        gene     = str(r['gene'])
        
        for idx, gt_type in enumerate(gt_types):
            if (gt_type == GT_HET or gt_type == GT_HOM_ALT):
                if gene != "None":
                    (key, value) = (idx_to_sample[idx], gene)
                    lof[idx_to_sample[idx]].append(gene)
    return lof
                
def sample_gene_interactions(c, args, idx_to_sample):
    #fetch variant gene dict for all samples
    get_variant_genes(c, args, idx_to_sample)
    #file handle for fetching the hprd graph
    file_graph = os.path.join(path_dirname, 'hprd_interaction_graph')
    #load the graph using cPickle
    gr = graph()
    gr = cPickle.load(open(file_graph, 'rb'))
    k = []
    variants = []
    #calculate nodes from the graph
    hprd_genes = gr.nodes()
    if args.gene == None or args.gene not in hprd_genes:
        sys.stderr.write("gene name either not provided or not represented \
                          in the p-p interaction file\n")
    elif args.gene in hprd_genes:
        x, y = breadth_first_search(gr,root=args.gene,filter=radius(args.radius))
        gst = digraph()
        gst.add_spanning_tree(x)
        dot = write(gst)
        out.write(dot)
        st, sd = shortest_path(gst, args.gene)
        
        if args.var_mode:
            for sample in sam.iterkeys():
                var = sam[str(sample)]
                #for each level return interacting genes if they are variants in the sample. 
                #0th order would be returned if the user chosen gene is a variant in the sample        
                for x in range(0, (args.radius+1)):
                    for each in var:
                        for key, value in sd.iteritems():
                            if value == x and key == each[0]:
                                print "\t".join([str(sample),str(args.gene), str(x), \
                                            str(key), str(each[1]), str(each[2]), str(each[3])])
                                            
        elif (not args.var_mode):
            for sample in sam.iterkeys():
                for each in sam[str(sample)]:
                    variants.append(each[0])
                for x in range(0, (args.radius+1)):
                    for key, value in sd.iteritems():
                        if value == x and key in set(variants):
                            k.append(key)
                    if k:
                        print "\t".join([str(sample), str(x)+"_order:", ",".join(k)])
                    else:
                        print "\t".join([str(sample), str(x)+"_order:", "none"])
                    #initialize keys for next iteration
                    k = []
                
            
def sample_lof_interactions(c, args, idx_to_sample):
    get_lof_genes(c, args, idx_to_sample)
    get_variant_genes(c, args, idx_to_sample)
    #file handle for fetching the hprd graph
    file_graph = os.path.join(path_dirname, 'hprd_interaction_graph')
    #load the graph using cPickle
    gr = graph()
    gr = cPickle.load(open(file_graph, 'rb'))
    #calculate nodes from the graph
    hprd_genes = gr.nodes()
    #initialize keys
    k = []
    variants = []
    
    if (not args.var_mode):
        for sample in lof.iterkeys():
            lofvariants = list(set(lof[str(sample)]))
            for each in sam[str(sample)]:
                variants.append(each[0])
                
            for gene in lofvariants:
                if gene in hprd_genes:
                    x, y = breadth_first_search(gr,root=gene,filter=radius(args.radius))
                    gst = digraph()
                    gst.add_spanning_tree(x)
                    st, sd = shortest_path(gst, gene)
                    #for each level return interacting genes if they are variants in the sample.
                    for rad in range(1, (args.radius+1)):
                        for key, value in sd.iteritems():
                            if (value == rad) and key in set(variants):
                                k.append(key)
                        if k:
                            print "\t".join([str(sample), str(gene), str(rad)+"_order:", ",".join(k)])
                        else:
                            print "\t".join([str(sample), str(gene), str(rad)+"_order:", "none"])
                        #initialize k
                        k = []
                else:
                    pass     
        
         
    elif args.var_mode:
        for sample in lof.iterkeys():
            lofvariants = list(set(lof[str(sample)]))
            var = sam[str(sample)]    
            for gene in lofvariants:
                if gene in hprd_genes:
                    x, y = breadth_first_search(gr,root=gene,filter=radius(args.radius))
                    gst = digraph()
                    gst.add_spanning_tree(x)
                    st, sd = shortest_path(gst, gene)
                    for rad in range(1, (args.radius+1)):
                        for each in var:
                            for key, value in sd.iteritems():
                                if value == rad and key == each[0]:
                                    print "\t".join([str(sample),str(gene), str(rad), str(key), \
                                                str(each[1]), str(each[2]), str(each[3])])
                else:
                    pass          
        
    
def sample_variants(c, args):
    idx_to_sample = util.map_indicies_to_samples(c)
    query = "SELECT DISTINCT variant_id, gt_types, gts, gene, impact, biotype \
             FROM variants"
    c.execute(query)
    print "\t".join(['sample','gene','order_of_interaction','interacting_gene', \
                     'var_id','impact','biotype'])
    sample_gene_interactions(c, args, idx_to_sample)
    
    
def sample_lof_variants(c, args):
    idx_to_sample = util.map_indicies_to_samples(c)
    query = "SELECT DISTINCT chrom, start, end, \
                             gt_types, gts, gene \
             FROM variants \
             WHERE is_lof='1'"
    c.execute(query)
      
    sample_lof_interactions(c, args, idx_to_sample) 
    
    
def variants(c, args):
    idx_to_sample = util.map_indicies_to_samples(c)
    query = "SELECT DISTINCT variant_id, gt_types, gts, gene, impact, biotype \
             FROM variants"
    c.execute(query)
    print "\t".join(['sample','lof_gene','order_of_interaction','interacting_gene', \
                     'var_id','impact','biotype'])
    get_variant_genes(c, args, idx_to_sample)
    
def genequery(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        sample_variants(c, args)

def lofgenequery(parser, args):
    if os.path.exists(args.db):
        conn = sqlite3.connect(args.db)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        variants(c, args)
        sample_lof_variants(c, args)