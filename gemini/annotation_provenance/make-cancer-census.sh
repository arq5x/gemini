DATE=`date "+%Y%m%d"`
wget -O - -q http://cancer.sanger.ac.uk/cancergenome/assets/cancer_gene_census.tsv \
| perl -pe 's/\r\n|\n|\r/\n/g' \
> cancer_gene_census.$DATE.tsv
