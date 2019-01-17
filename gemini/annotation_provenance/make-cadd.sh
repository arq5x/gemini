# Download the raw CADD TSV and Tabix index (no annotations, just scores)
<<DOCS
wget http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz.tbi

# E.g snippet of contents in the above file
ChromPosRefAltRawScorePHRED
110002AC-0.6412591.170
17010002AG-0.6473071.149
14910002AT-0.5333261.574

#compress columns 5 & 6 to collapse values for all possible 'Alt' for a given position
# as (lexically ordered) comma separated values, using Bedtools
# compress float values to 2 decimal precision

#resultant compression would look like this
ChromPosRefAltRawScorePHREDPosRefRawScorePHRED
14910002AT10002A-0.64,-0.65,-0.53332611.17,1.15,1.57


DOCS
set -euo pipefail

#wget -c https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz 
wget -c http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/InDels.tsv.gz

zgrep -v ^# whole_genome_SNVs.tsv.gz \
	    | awk '{printf("%s\t%d\t%s\t%s\t%0.2f\t%0.2f\n",$1,$2,$3,$4,$5,$6)}' \
	    | bedtools groupby -g 1,2,3 -c 5,6 -o collapse,collapse \
		| sort -k1,1 -k2,2n \
		| bgzip -l 7 -c > caddv1.4.compressed.gz
tabix -b 2 -e 2 caddv1.4.compressed.gz
