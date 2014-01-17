####################################################
# For a detailed gene table and a summary gene table
####################################################


#!/usr/bin/env python

import sys
import os
import itertools
from collections import defaultdict
from sets import Set


filename = 'detailed_gene_table'
detailed_out = open(filename, 'w')

file = 'summary_gene_table'
summary_out = open(file, 'w')

# write out files for detailed and summary gene table
detailed_out.write("\t".join(["Chromosome","Gene_name","Is_hgnc","Ensembl_gene_id","Ensembl_transcript_id","Biotype","Transcript_status",
                     "CCDS_id","HGNC_id","CDS_length","Protein_length","Transcript_start","Transcript_end","strand","Synonyms", "Rvis_pct"]))     
detailed_out.write("\n")

summary_out.write("\t".join(["Chromosome","Gene_name","Is_hgnc","Ensembl_gene_id",
                         "HGNC_id","Synonyms", "Rvis_pct","Strand","Transcript_min_start","Transcript_max_end"]))     
summary_out.write("\n")



genic_intolerance = defaultdict(list)
keygene = list_hgnc = []

#initializing values for the summary gene table
transcript_min = defaultdict(list)
transcript_max = defaultdict(list)
lines_seen = set()


for line in open("genic_intolerance_dataset2", 'r'):
    if line.startswith("#") is False:
        field = line.strip().split("\t")
        name = str(field[0])
        score =  str(field[1])
        percentile = str(field[2])
        (key,value) = (name, percentile)
        genic_intolerance[name].append(percentile)
        

# Dictionary for summary gene table to handle transcript min, max co-ordinates
for each in open("raw_gene_table", 'r'):
    if each.startswith("Chromosome") is False:
        k = each.strip().split("\t")
        chr = "chr"+str((k[0]))
        ens = str(k[2])
        start = str(k[10])
        end = str(k[11])
        transcript_min[(chr,ens)].append(start)
        transcript_max[(chr,ens)].append(end)



for each in open("raw_gene_table", 'r'):
    if each.startswith("Chromosome") is False:
        k = each.strip().split("\t")
        chrom = "chr"+str((k[0]))
        hgnc = str(k[1])
        ens_geneid = str(k[2])
        ens_transid = str(k[3])
        trans_biotype = str(k[4])
        status = str(k[5])
        ccds_id = str(k[6]) #these id's are unique to transcripts
        hgnc_id = str(k[7])
        cds_len = str(k[8])
        protein_len = str(k[9])
        transcript_start = str(k[10])
        transcript_end = str(k[11])
        strand = str(k[12])
        #remove space between names
        previous =  str(k[13]).replace(" ","")
        synonyms = str(k[14]).replace(" ","")
        
        # sort all transcript start and end positions for a gene (use ens_geneid, since HGNC is not always true)
        # Capture the first and the last position from the sorted list to give min, max
        if (chrom,ens_geneid) in transcript_min:
            minmum = sorted(transcript_min[(chrom,ens_geneid)])[0]
        if (chrom,ens_geneid) in transcript_max:
            maxmum = sorted(transcript_max[(chrom,ens_geneid)])[-1]
        
        rvis = genic_intolerance[hgnc][0] if hgnc in genic_intolerance else None
        
        if hgnc != "None":
            list_hgnc.append(hgnc)
        #we don't want string of Nones
        if "None" in previous and "None" in synonyms and "None" in hgnc:
            string = None
        else:
            # We would like all genes names to be put together
            gene_string = hgnc+","+previous+","+synonyms
            
            #get rid of Nones in gene strings
            if gene_string.startswith("None"):
                string = gene_string.replace("None,","")
            else:
                string = gene_string.replace(",None","")

        #Nonetype object has no attribute split
        if string is not None:
            genes = set(string.split(","))
            if len(genes) > 1:
            # We would like to represent each member of the gene list as a key and the remainder as synonyms each time
                for each in genes:
                    keygene = set([each])
                    synonym = genes.difference(keygene)
                    gene_name = ','.join(keygene)
                    other_names = ','.join(synonym)
                    hgnc_flag = "1" if gene_name in list_hgnc else "0"
                    # only when the gene is a HGNC name, it would have an hgnc id
                    is_hgnc_id = hgnc_id if gene_name in list_hgnc else "None"
                    
                    # handling duplicate lines (due to transcripts) in summary table (which we don't care for in this table)
                    # writing to outfile for the summary gene table
                    line = "\t".join([chrom,gene_name,hgnc_flag,ens_geneid,is_hgnc_id,
                                         other_names,str(rvis),strand,minmum,maxmum])
                    if line not in lines_seen:
                        summary_out.write(line)
                        summary_out.write("\n")
                        lines_seen.add(line)
                    
                    
                    # Writing to out for detailed gene table
                    detailed_out.write("\t".join([chrom,gene_name,hgnc_flag,ens_geneid,ens_transid,trans_biotype,status,ccds_id,is_hgnc_id,cds_len,
                                         protein_len,transcript_start,transcript_end,strand,other_names,str(rvis)]))
                    detailed_out.write("\n")
            
            # if there is one gene name in the list, we just want it to be the key
            elif len(genes) == 1:
                gene_name = ','.join(genes)
                other_names = "None"
                hgnc_flag = "1" if gene_name in list_hgnc else "0"
                is_hgnc_id = hgnc_id if gene_name in list_hgnc else "None"
                
                
                # handling duplicate lines (due to transcripts) in summary table (which we don't care for in this table)
                # writing to outfile for the summary gene table
                line = "\t".join([chrom,str(gene_name),hgnc_flag,ens_geneid,is_hgnc_id,
                                  other_names,str(rvis),strand,minmum,maxmum])
                                  
                if line not in lines_seen:
                    summary_out.write(line)
                    summary_out.write("\n")
                    lines_seen.add(line)
                
                
                # write to out for detailed gene table
                detailed_out.write("\t".join([chrom,str(gene_name),hgnc_flag,ens_geneid,ens_transid,trans_biotype,status,ccds_id,is_hgnc_id,cds_len,
                                               protein_len,transcript_start,transcript_end,strand,other_names,str(rvis)]))
                detailed_out.write("\n")
        # if there are no HGNC, previous or synonyms names for an ensembl entry, just return None
        elif string is None:
            gene_name = "None"
            other_names = "None"
            hgnc_flag = "0"
            is_hgnc_id = "None"
            
            #handling duplicate lines (due to transcripts) in summary table (which we don't care for in this table)
            #writing to outfile for the summary gene table
            line = "\t".join([chrom,gene_name,hgnc_flag,ens_geneid,is_hgnc_id,
                              other_names,str(rvis),strand,minmum,maxmum])
            if line not in lines_seen:
                summary_out.write(line)
                summary_out.write("\n")
                lines_seen.add(line)
            
            # probably we still want to print these lines where gene is none since ensembl gene id has value
            detailed_out.write("\t".join([chrom,gene_name,hgnc_flag,ens_geneid,ens_transid,trans_biotype,status,ccds_id,is_hgnc_id,cds_len,
                                           protein_len,transcript_start,transcript_end,strand,other_names,str(rvis)]))
            detailed_out.write("\n")
            
detailed_out.close()
summary_out.close()
 