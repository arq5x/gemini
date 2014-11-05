import sys
import os
import database
import database_postgresql
from gemini.config import read_gemini_config

class gene_detailed:

    def __init__(self, fields):
        
        self.fields = [f if f != "None" else None for f in fields]
        self.chrom = self.fields[0]
        self.gene = self.fields[1]
        self.is_hgnc = self.fields[2]
        self.ensembl_gene_id = self.fields[3]
        self.ensembl_trans_id = self.fields[4]
        self.biotype = self.fields[5]
        self.trans_status = self.fields[6]
        self.ccds_id = self.fields[7]
        self.hgnc_id = self.fields[8]
        self.cds_length = self.fields[9]
        self.protein_length = self.fields[10]
        self.transcript_start = self.fields[11]
        self.transcript_end = self.fields[12]
        self.strand = self.fields[13]
        self.synonym = self.fields[14]
        self.rvis = self.fields[15]
        self.entrez = self.fields[16]
        self.mam_phenotype = self.fields[17]
        
    def __str__(self):
        return ",".join([self.chrom, self.gene, self.is_hgnc, self.ensembl_gene_id, self.ensembl_trans_id, self.biotype, self.trans_status,
                         self.ccds_id, self.hgnc_id, self.cds_length, self.protein_length, self.transcript_start, self.transcript_end, self.strand, 
                         str(self.synonym), self.rvis, self.entrez, self.mam_phenotype])
                                             
class gene_summary:
    
    def __init__(self, fields):
         self.fields = [f if f != "None" else None for f in fields]
         self.columns = self.fields[:]
         self.chrom = self.fields[0]
         self.gene = self.fields[1]
         self.is_hgnc = self.fields[2]
         self.ensembl_gene_id = self.fields[3]
         self.hgnc_id = self.fields[4]
         self.synonym = self.fields[5]
         self.rvis = self.fields[6]
         self.strand = self.fields[7]
         self.transcript_min_start = self.fields[8]
         self.transcript_max_end = self.fields[9]
         self.mam_phenotype = self.fields[10]
         
    def __str__(self):
        return ",".join([self.chrom, self.gene, self.is_hgnc, self.ensembl_gene_id, self.hgnc_id, self.synonym, self.rvis, 
                         self.strand, self.transcript_min_start, self.transcript_max_end, self.mam_phenotype])
         
def update_cosmic_census_genes( cursor, args ):
    """
    Update the gene summary table with
    whether or not a given gene is in the
    COSMIC cancer gene census
    """
    config = read_gemini_config( args= args )
    path_dirname = config["annotation_dir"]
    file = os.path.join(path_dirname, 'cancer_gene_census.20140120.tsv')
    
    cosmic_census_genes = []
    for line in open(file, 'r'):
        fields = line.strip().split("\t")
        gene = fields[0]
        chrom = "chr" + fields[3]
        cosmic_census_genes.append((1,gene,chrom))

    if args.dbtype == "sqlite":
        database.update_gene_summary_w_cancer_census(cursor, cosmic_census_genes)
    elif args.dbtype == "postgresql":
        database_postgresql.update_gene_summary_w_cancer_census(cursor, cosmic_census_genes)
         
        
