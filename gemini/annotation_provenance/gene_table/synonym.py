#!/usr/bin/env python

import sys
import os
import itertools
from collections import defaultdict

files = 'gene_table'
outfile = open(files, 'w')

outfile.write("\t".join(["Chromosome","HGNC_symbol","Ensembl_gene_id","Ensembl_transcript_id","Biotype","Transcript_status",
                         "CCDS_id","HGNC_id","CDS_length","Protein_length","transcript_start","transcript_end","strand","Previous_symbol","Synonymous"]))
outfile.write("\n")

hgncfile = defaultdict(list)

for lines in open("hgnc_file",'r'):
    if lines.startswith("HGNC") is False:
        seq = lines.strip().split("\t")
        (key,value) = (seq[0],(seq[2],seq[3]))
        hgncfile[key].append(value)
    
with open("ensembl_format") as f:
    for each in f:
        if each.startswith("Chromosome") is False:
            field = each.strip().split("\t")
            hgncid = field[7]
            if hgncid in hgncfile:
                for values in hgncfile[hgncid]:
                    #print values[0], values[1]
                    previous_name = values[0]
                    synonymous = values[1]
                sequence = [field[0],field[1],field[2],field[3],field[4],field[5],field[6],field[7],field[8],field[9],field[10],field[11],field[12],previous_name, synonymous]
            else:
                # withdrawn entries (e.g 40184) or genes that do not have an HGNC symbol (None) will not have entry in hgnc_file
                previous_name = "None"
                synonymous = "None"
                sequence = [field[0],field[1],field[2],field[3],field[4],field[5],field[6],field[7],field[8],field[9],field[10],field[11],field[12],previous_name, synonymous]
            
            outfile.write("\t".join(sequence))
            outfile.write("\n")                
            
outfile.close()
