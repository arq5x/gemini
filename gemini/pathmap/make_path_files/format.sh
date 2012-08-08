# Pathways mapped for ensembl genes based on gene names and uniprot ids

############################################################################################################################################
# Processing the ensembl external id map file from biomart (with transcript column) mart_export_transcripts.txt (creation date:may31st 2012)
# build ensembl genes 66
############################################################################################################################################
cat mart_export_transcripts66.txt | awk 'NR==1' > mart_header
cat mart_export_transcripts66.txt | awk 'NR==2, EOF'| cut -f1,2 | sed '/^$/d' | awk '{print$1}' > file1.1 #this basically merges the two columns(swissprot/trembl) such that swissprot is chosen when both ids are present
cat mart_export_transcripts66.txt | awk 'NR==2, EOF' | cut -f3,4,5,6 > file2.2
paste file1.1 file2.2 | awk -F'\t' '{ OFS = "\t" }; {for(n=1; n<=NF; n++) sub(/^$/, "none", $n); print $0}' > ensembl_genes
cat ensembl_genes | sort | uniq > ensembl_genes66
rm file1.1 file2.2 ensembl_genes
###########################################################
# Formatting kegg pathcode/pathname file
###########################################################
cat KEGG_REST_files/pathways | cut -f1 > kegg_path
cat KEGG_REST_files/pathways | cut -f 2- | sed 's/ - Homo sapiens (human)//g' | sed 's/ \/ /_/g' | sed 's/ (/_/g' | sed 's/)//g' | sed 's/, /_/g' | sed 's/ - /_/g' | sed 's/-/_/g' | sed 's/ /_/g' > pathways-fmt
paste  kegg_path pathways-fmt > kegg_pathway
rm kegg_path pathways-fmt

################################################
# Formatting the hsa to uniprot map file
################################################
cat KEGG_REST_files/all_kegg_hsa_to_uniprot | sed 's/^up://g' > uniprot_hsa_map

######################################
# Processing the (hsa to gene name) file
########################################

cat KEGG_REST_files/kegg_hsa_to_gene | cut -d";" -f1 > kegg_hsa_to_gene_fmt

########################################
# python code to generate the mypathfile
########################################

python unihsagene.py

Establishes a one to one relationship across;
uniprot, gene (agn,hgnc,ensemblid,transcript), hsa, path



