# Download the raw CADD TSV and Tabix index (no annotations, just scores)
wget http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.0/whole_genome_SNVs.tsv.gz.tbi

# E.g snippet of contents in the above file
Chrom	Pos	Ref	Alt	RawScore	PHRED
1	10002	A	C	-0.641259	1.170
1	10002	A	G	-0.647307	1.149
1	10002	A	T	-0.533326	1.574

#compress columns 5 & 6 to collapse values for all possible 'Alt' for a given position
# as (lexically ordered) comma separated values, using Bedtools
# compress float values to 2 decimal precision

#resultant compression would look like this
Chrom	Pos	Ref	RawScore	PHRED
1	10002	A	-0.64,-0.65,-0.53	1.17,1.15,1.57


#compression
zcat whole_genome_SNVs.tsv.gz | grep -v "^#" | awk '{printf("%s\t%d\t%s\t%s\t%0.2f\t%0.2f\n",$1,$2,$3,$4,$5,$6)}' \
| bedtools groupby -g 1,2,3 -c 5,6 -o collapse,collapse | bgzip > whole_genome_SNVs.tsv.compressed.gz

#index
tabix -b 2 -e 2 whole_genome_SNVs.tsv.compressed.gz

