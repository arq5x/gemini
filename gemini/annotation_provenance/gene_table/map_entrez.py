#!/usr/bin/env python

import sys
import os
import itertools
from collections import defaultdict

files = 'raw_gene_table'
outfile = open(files, 'w')

outfile.write("\t".join(["Chromosome","HGNC_symbol","Ensembl_gene_id","Ensembl_transcript_id","Biotype","Transcript_status",
                         "CCDS_id","HGNC_id","CDS_length","Protein_length","transcript_start","transcript_end",
                         "strand","Previous_symbol","Synonymous","entrez_id"]))
outfile.write("\n")

entrez = defaultdict(list)

for lines in open("ensembl75_3",'r'):
    if lines.startswith("Ensembl") is False:
        seq = lines.strip().split("\t")
        (key,value) = (seq[1],seq[2])
        entrez[key].append(value)
  
with open("gene_table") as f:
    for each in f:
        if each.startswith("Chromosome") is False:
            field = each.strip().split("\t")
            transcript = field[3]
            if transcript in entrez:
                for value in entrez[transcript]:
                    sequence = [field[0],field[1],field[2],field[3],field[4],field[5],field[6],field[7],field[8],field[9],field[10],field[11],field[12],field[13],field[14],value]
            else:
                # return none for entrez where there is no mapping
                value = "None"
                sequence = [field[0],field[1],field[2],field[3],field[4],field[5],field[6],field[7],field[8],field[9],field[10],field[11],field[12],field[13],field[14],value]
            
            outfile.write("\t".join(sequence))
            outfile.write("\n")                
            
outfile.close()
