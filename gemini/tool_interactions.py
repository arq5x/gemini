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

def get_variant_genes(c, args, idx_to_sample):
    samples = defaultdict(list)
    for r in c:
        gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
        gts      = np.array(cPickle.loads(zlib.decompress(r['gts'])))
        var_id = str(r['variant_id'])
        chrom = str(r['chrom'])
        start = str(r['start'])
        end = str(r['end'])
        gene     = str(r['gene'])
        impact = str(r['impact'])
        biotype = str(r['biotype'])
        in_dbsnp = str(r['in_dbsnp'])
        clinvar_sig = str(r['clinvar_sig'])
        clinvar_disease_name = str(r['clinvar_disease_name'])
        aaf_1kg_all = str(r['aaf_1kg_all'])
        aaf_esp_all = str(r['aaf_esp_all'])

        for idx, gt_type in enumerate(gt_types):
            if (gt_type == HET or gt_type == HOM_ALT):
                if gene != "None":
                    (key, value) = (idx_to_sample[idx], \
                                   (gene,var_id,chrom,start,end,impact, \
                                   biotype,in_dbsnp,clinvar_sig, \
                                   clinvar_disease_name,aaf_1kg_all, \
                                   aaf_esp_all))
                    samples[idx_to_sample[idx]].append(value)
    return samples

def get_lof_genes(c, args, idx_to_sample):
    lof = defaultdict(list)
    for r in c:
        gt_types = np.array(cPickle.loads(zlib.decompress(r['gt_types'])))
        gts      = np.array(cPickle.loads(zlib.decompress(r['gts'])))
        gene     = str(r['gene'])

        for idx, gt_type in enumerate(gt_types):
            if (gt_type == HET or gt_type == HOM_ALT):
                if gene != "None":
                    (key, value) = (idx_to_sample[idx], gene)
                    lof[idx_to_sample[idx]].append(gene)
    return lof

def sample_gene_interactions(c, args, idx_to_sample):
    out = open("file.dot", 'w')
    #fetch variant gene dict for all samples
    samples = get_variant_genes(c, args, idx_to_sample)
    #file handle for fetching the hprd graph
    config = read_gemini_config( args = args )
    path_dirname = config["annotation_dir"]
    file_graph = os.path.join(path_dirname, 'hprd_interaction_graph')
    #load the graph using cPickle and close file handle
    gr = graph()
    f = open(file_graph, 'rb')
    gr = cPickle.load(f)
    f.close()
    k = []
    variants = []
    #calculate nodes from the graph
    hprd_genes = gr.nodes()
    if args.gene == None or args.gene not in hprd_genes:
        sys.stderr.write("Gene name not found or")
        sys.stderr.write(" gene not in p-p interaction file\n")

    elif args.gene in hprd_genes:
        x, y = \
            breadth_first_search(gr,root=args.gene,filter=radius(args.radius))
        gst = digraph()
        gst.add_spanning_tree(x)
        dot = write(gst)
        out.write(dot)
        st, sd = shortest_path(gst, args.gene)

        if args.var_mode:
            for sample in samples.iterkeys():
                var = samples[str(sample)]
                #for each level return interacting genes if they are
                # variants in the sample.
                # 0th order would be returned if the user chosen
                # gene is a variant in the sample
                for x in range(0, (args.radius+1)):
                    for each in var:
                        for key, value in sd.iteritems():
                            if value == x and key == each[0]:
                                print "\t".join([str(sample),str(args.gene), \
                                          str(x), \
                                          str(key), \
                                          str(each[1]), \
                                          str(each[2]), \
                                          str(each[3]), \
                                          str(each[4]), \
                                          str(each[5]), \
                                          str(each[6]), \
                                          str(each[7]), \
                                          str(each[8]), \
                                          str(each[9]), \
                                          str(each[10]), \
                                          str(each[11])])
        elif (not args.var_mode):
            for sample in samples.iterkeys():
                for each in samples[str(sample)]:
                    variants.append(each[0])
                for x in range(0, (args.radius+1)):
                    for key, value in sd.iteritems():
                        if value == x and key in set(variants):
                            k.append(key)
                    if k:
                        print "\t".join([str(sample), str(args.gene), \
                                 str(x)+"_order:",
                                 ",".join(k)])
                    else:
                        print "\t".join([str(sample), str(args.gene), \
                                         str(x)+"_order:", "none"])
                    #initialize keys for next iteration
                    k = []
                #initialize variants list for next iteration
                variants = []


def sample_lof_interactions(c, args, idx_to_sample, samples):
    lof = get_lof_genes(c, args, idx_to_sample)
    #file handle for fetching the hprd graph
    config = read_gemini_config( args = args )
    path_dirname = config["annotation_dir"]
    file_graph = os.path.join(path_dirname, 'hprd_interaction_graph')
    #load the graph using cPickle and close file handle
    gr = graph()
    f = open(file_graph, 'rb')
    gr = cPickle.load(f)
    f.close()
    #calculate nodes from the graph
    hprd_genes = gr.nodes()
    #initialize keys
    k = []
    variants = []

    if (not args.var_mode):
        for sample in lof.iterkeys():
            lofvariants = list(set(lof[str(sample)]))
            for each in samples[str(sample)]:
                variants.append(each[0])
            for gene in lofvariants:
                if gene in hprd_genes:
                    x, y = \
                        breadth_first_search(gr,root=gene,\
                        filter=radius(args.radius))

                    gst = digraph()
                    gst.add_spanning_tree(x)
                    st, sd = shortest_path(gst, gene)
                    # for each level return interacting genes
                    # if they are variants in the sample.
                    for rad in range(1, (args.radius+1)):
                        for key, value in sd.iteritems():
                            if (value == rad) and key in set(variants):
                                k.append(key)
                        if k:
                            print "\t".join([str(sample), \
                                       str(gene), \
                                       str(rad)+"_order:",
                                       ",".join(k)])
                        else:
                            print "\t".join([str(sample), \
                                       str(gene), \
                                       str(rad)+"_order:", \
                                       "none"])
                        #initialize k
                        k = []
            #initialize variants list for next iteration
            variants = []
    elif args.var_mode:
        for sample in lof.iterkeys():
            lofvariants = list(set(lof[str(sample)]))
            var = samples[str(sample)]
            for gene in lofvariants:
                if gene in hprd_genes:
                    x, y = \
                         breadth_first_search(gr,root=gene, \
                         filter=radius(args.radius))
                    gst = digraph()
                    gst.add_spanning_tree(x)
                    st, sd = shortest_path(gst, gene)
                    for rad in range(1, (args.radius+1)):
                        for each in var:
                            for key, value in sd.iteritems():
                                if value == rad and key == each[0]:
                                    print "\t".join([str(sample), \
                                               str(gene), \
                                               str(rad), \
                                               str(key), \
                                               str(each[1]), \
                                               str(each[2]), \
                                               str(each[3]), \
                                               str(each[4]), \
                                               str(each[5]), \
                                               str(each[6]), \
                                               str(each[7]), \
                                               str(each[8]), \
                                               str(each[9]), \
                                               str(each[10]), \
                                               str(each[11])])


def sample_variants(c, args):
    idx_to_sample = util.map_indices_to_samples(c)
    query = "SELECT variant_id, gt_types, gts, gene, impact, biotype, \
                    in_dbsnp, clinvar_sig, clinvar_disease_name, aaf_1kg_all, aaf_esp_all, chrom, \
                    start, end  \
             FROM variants"
    c.execute(query)

    if args.command == 'interactions':
        #header
        if args.var_mode:
            print "\t".join(['sample','gene','order_of_interaction', \
                             'interacting_gene', 'var_id', 'chrom', 'start', \
                             'end', 'impact', 'biotype', 'in_dbsnp', \
                             'clinvar_sig', 'clinvar_disease_name', 'aaf_1kg_all', \
                             'aaf_esp_all'])

        if (not args.var_mode):
            print "\t".join(['sample','gene','order_of_interaction', \
                     'interacting_gene'])
        sample_gene_interactions(c, args, idx_to_sample)

    elif args.command == 'lof_interactions':
        samples = get_variant_genes(c, args, idx_to_sample)
        return samples


def sample_lof_variants(c, args, samples):
    idx_to_sample = util.map_indices_to_samples(c)
    query = "SELECT chrom, start, end, \
                             gt_types, gts, gene \
             FROM variants \
             WHERE is_lof='1'"
    c.execute(query)

    #header
    if args.var_mode:
        print "\t".join(['sample','lof_gene','order_of_interaction', \
                    'interacting_gene', 'var_id', 'chrom', 'start', \
                    'end', 'impact','biotype','in_dbsnp', 'clinvar_sig', \
                    'clinvar_disease_name', 'aaf_1kg_all','aaf_esp_all'])

    elif (not args.var_mode):
        print "\t".join(['sample','lof_gene','order_of_interaction', \
                         'interacting_gene'])

    sample_lof_interactions(c, args, idx_to_sample, samples)


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
        samples = sample_variants(c, args)
        sample_lof_variants(c, args, samples)
