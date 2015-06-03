#!/usr/bin/env python

from collections import defaultdict

filename = 'ensembl_format'
out = open(filename, 'w')

out.write("\t".join(["Chromosome","HGNC_symbol","Ensembl_gene_id","Ensembl_transcript_id","Biotype","Transcript_status",
                     "CCDS_id","HGNC_id","CDS_length","Protein_length","transcript_start","transcript_end","strand"]))     
out.write("\n")

ccdslen = defaultdict(list)
ensemb = defaultdict(list)


for each in open("ensembl75_2", 'r'):
    if each.startswith("Ensembl") is False:
        col = each.strip().split("\t")
        if col[2] != "None":
            Protein_length = int(col[2])/int(3) - 1
        else:
            Protein_length = "None"
        (key, value) = (col[1], (col[2],str(Protein_length)))
        ccdslen[key].append(value)

for line in open("ensembl75_1", 'r'):
    if line.startswith("Chromosome") is False:
        k = line.strip().split("\t")
        chrom = str((k[0]))
        hgnc = str(k[1])
        ens_geneid = str(k[2])
        ens_transid = str(k[3])
        trans_biotype = str(k[4])
        status = str(k[5])
        ccds_id = str(k[6])
        hgnc_id = str(k[7])
        transcript_start = str(k[8])
        transcript_end = str(k[9])
        strand = str(k[10])
    
        if ens_transid in ccdslen:
            for each in ccdslen[ens_transid]:
                cds_len = each[0]
                protein_len = each[1]
            string = [chrom,hgnc,ens_geneid,ens_transid,trans_biotype,status,ccds_id,hgnc_id,cds_len,protein_len,transcript_start,transcript_end,strand]
        else:
            print "line fail"
        out.write("\t".join(string))
        out.write("\n")
        
out.close()
