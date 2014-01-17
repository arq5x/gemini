#!/usr/bin/env python

import sys

class gene_detailed:

    def __init__(self, field):
        
        self.fields = field[:]
        self.chrom = field[0]
        self.gene = field[1]
        self.is_hgnc = field[2]
        self.ensembl_gene_id = field[3]
        self.ensembl_trans_id = field[4]
        self.biotype = field[5]
        self.trans_status = field[6]
        self.ccds_id = field[7]
        self.hgnc_id = field[8]
        self.cds_length = field[9]
        self.protein_length = field[10]
        self.transcript_start = field[11]
        self.transcript_end = field[12]
        self.strand = field[13]
        self.synonym = field[14]
        self.rvis = field[15]
        
        
        
    def __str__(self):
        return ",".join([self.chrom, self.gene, self.is_hgnc, self.ensembl_gene_id, self.ensembl_trans_id, self.biotype, self.trans_status,
                         self.ccds_id, self.hgnc_id, self.cds_length, self.protein_length, self.transcript_start, self.transcript_end, self.strand, 
                         str(self.synonym), self.rvis])
                     
                     
class gene_summary:
    
    def __init__(self, col):
         self.columns = col[:]
         self.chrom = col[0]
         self.gene = col[1]
         self.is_hgnc = col[2]
         self.ensembl_gene_id = col[3]
         self.hgnc_id = col[4]
         self.synonym = col[5]
         self.rvis = col[6]
         self.strand = col[7]
         self.transcript_min_start = col[8]
         self.transcript_max_end = col[9]
         
    def __str__(self):
        return ",".join([self.chrom, self.gene, self.is_hgnc, self.ensembl_gene_id, self.hgnc_id, self.synonym, self.rvis, 
                         self.strand, self.transcript_min_start, self.transcript_max_end])
         
         
          
         
        
