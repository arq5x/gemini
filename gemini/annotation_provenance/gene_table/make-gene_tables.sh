#!/bin/bash

set -euo pipefail

# get biomart perl API
git clone --branch cvs/release-0_7 https://github.com/biomart/biomart-perl

# get Ensembl biomart registry
wget -O ensembl_biomart.xml https://www.ensembl.org/biomart/martservice?type=registry

# get Ensembl version from the downloaded registry xml
ENSEMBLVERSION=`cat ensembl_biomart.xml | grep ensembl_mart | perl -pe 's/.*ensembl_mart_(\d+).*/$1/'`

# initialize biomart cache - otherwise messages about unavailable cache go to
# stdout and mess up the tables below
perl -Ibiomart-perl/lib downloads/init_biomart.pl

# get Ensembl gene tables
perl -Ibiomart-perl/lib downloads/query1.pl > mart_export1
perl -Ibiomart-perl/lib downloads/query2.pl > mart_export2
perl -Ibiomart-perl/lib downloads/query3.pl > mart_export3

# get HGNC gene symbols
perl downloads/hgnc.query.pl > HGNC_download

# post-process downloaded files
cat mart_export1 | awk -F'\t' '{ OFS = "\t" }; {for(n=1; n<=NF; n++) sub(/^$/, "None", $n); print $0}' > ensembl_1
cat mart_export2 | awk -F'\t' '{ OFS = "\t" }; {for(n=1; n<=NF; n++) sub(/^$/, "None", $n); print $0}' > ensembl_2
cat mart_export3 | awk -F'\t' '{ OFS = "\t" }; {for(n=1; n<=NF; n++) sub(/^$/, "None", $n); print $0}' > ensembl_3
cat HGNC_download | awk -F'\t' '{ OFS = "\t" }; {for(n=1; n<=NF; n++) sub(/^$/, "None", $n); print $0}' > hgnc_file
cat anno_files/HMD_HumanPhenotype.rpt | awk -F'\t' '{ OFS = "\t" }; {for(n=1; n<=NF; n++) sub(/^$/, "None", $n); print $0}' > HMD_HumanPhenotype

# make ensembl_format table
python ensembl.py

# make gene_table
python synonym.py

# make raw_gene_table
python map_entrez.py

# make tables detailed_gene_table_v$ENSEMBLVERSION and summary_gene_table_v$ENSEMBLVERSION
python combined_gene_table.py
mv detailed_gene_table detailed_gene_table_v$ENSEMBLVERSION
mv summary_gene_table summary_gene_table_v$ENSEMBLVERSION

# remove temp files and cloned biomart-perl repo and caches
rm gene_table raw_gene_table ensembl_* hgnc_file HMD_HumanPhenotype HGNC_download mart_export*
rm -rf biomart-perl
rm -rf Cached cachedRegistries

