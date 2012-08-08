#!/usr/bin/env python

import os
import sys
import sqlite3
from collections import defaultdict

path_dirname = '/Users/udp3f/python_prog/gemini/gemini/pathmap'

def get_pathfile(args):     
        
    if args.version == None or args.version == 'v66':
        path_file = os.path.join(path_dirname, 'mypathfile66') # default path map file set to 66
    elif args.version == 'v67':
        path_file = os.path.join(path_dirname, 'mypathfile67')
    elif args.version == 'v68':
        path_file = os.path.join(path_dirname, 'mypathfile68')
    else:
        sys.stderr.write("Unsupported version\n")
    
    if args.version in (None,'v66', 'v67', 'v68'):       
        pathway = defaultdict(list)
        for line in open(path_file, 'r'):
            fields=line.strip().split("\t")
            uniprot = fields[0]
            agn = fields[1]
            hgnc = fields[2]
            ensid = fields[3]
            ensript = fields[4]
            hsa = fields[5]
            path = fields[6]
            (key, value) = ((agn,hgnc,ensid,ensript), path)
            pathway[key].append(value) 
        return pathway


def get_all_pathways(c, args):
    
    pathway = get_pathfile(args)
    if pathway:
        query = "SELECT DISTINCT gene, transcript from variant_impacts"
        c.execute(query)
        for row in c:
            gene = str(row[0])
            transcript = str(row[1])
    
            p = [(gene, transcript, v) for k, v in pathway.iteritems() if (k[0] == gene and k[3] == transcript) or (k[1] == gene and k[3] == transcript) or (k[2] == gene and k[3] == transcript)]
            if p:
                count = 0
                for each in p:
                    for path in each[2]:  
                        if path == "None":
                            count +=1
                    if len(each[2]) == count:
                        print each[0], each[1], None
                    elif count == 0:
                        print each[0], each[1], ','.join(each[2])
                    if count >= 1 and len(each[2]) != count:
                        for items in each[2]:
                            if items != "None":
                                print each[0], each[1], items
            else:
                print gene, transcript, None  #None for genes that are not represented in the kegg pathfile e.g MIR34A


def get_lof_pathways(c, args):
    
    pathway = get_pathfile(args)
    if pathway:
        query = "SELECT DISTINCT gene, transcript from variant_impacts where is_lof='1'"
        c.execute(query)
        for each in c:
            af_gene = str(each[0])
            af_transcript = str(each[1])
        
            p = [(af_gene, af_transcript, v) for k, v in pathway.iteritems() if (k[0] == af_gene and k[3] == af_transcript) or (k[1] == af_gene and k[3] == af_transcript) or (k[2] == af_gene and k[3] == af_transcript)]
            if p:
                count = 0
                for each in p:
                    for path in each[2]:  
                        if path == "None":
                            count +=1
                    if len(each[2]) == count:
                        #print each[0], each[1], None
                        pass
                    elif count == 0:
                        print each[0], each[1], ','.join(each[2])
                    if count >= 1 and len(each[2]) != count:
                        for items in each[2]:
                            if items != "None":
                                print each[0], each[1], items
            #else:
            #    print af_gene, af_transcript, None  #None for genes that are not represented in the kegg pathfile e.g MIR34A        


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
        
        
                    
                    
                                