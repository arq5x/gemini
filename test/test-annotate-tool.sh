###########################################################################################
#1. Test annotating variants using the "boolean" function
###########################################################################################
echo "    annotate-tool.t1...\c"

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548
chr1	30922	30923" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -c anno -t boolean test.snpeff.vcf.db

echo "chr1	30548	1
chr1	30860	0
chr1	30869	0
chr1	30895	0
chr1	30923	1
chr1	69270	0
chr1	69428	0
chr1	69511	0
chr1	69761	0
chr1	69871	0" > exp

gemini query -q "select chrom, end, anno from variants" \
	test.snpeff.vcf.db > obs
check obs exp
rm obs exp
rm *.gz*


###########################################################################################
#2. Test annotating variants using the "count" function
###########################################################################################
echo "    annotate-tool.t2...\c"

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548
chr1	30920	30925
chr1	30922	30923" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -c anno2 -t count test.snpeff.vcf.db

echo "chr1	30548	1
chr1	30860	0
chr1	30869	0
chr1	30895	0
chr1	30923	2
chr1	69270	0
chr1	69428	0
chr1	69511	0
chr1	69761	0
chr1	69871	0" > exp

gemini query -q "select chrom, end, anno2 from variants" \
	test.snpeff.vcf.db > obs
check obs exp
rm obs exp
rm *.gz*


###########################################################################################
#3. Test annotating variants using the "list" function
###########################################################################################
echo "    annotate-tool.t3...\c"

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a
chr1	30920	30925	b
chr1	30922	30923	c" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -c anno3 -t list -e 4 test.snpeff.vcf.db

echo "chr1	30548	a
chr1	30860	None
chr1	30869	None
chr1	30895	None
chr1	30923	b,c
chr1	69270	None
chr1	69428	None
chr1	69511	None
chr1	69761	None
chr1	69871	None" > exp

gemini query -q "select chrom, end, anno3 from variants" \
	test.snpeff.vcf.db > obs
check obs exp
rm obs exp
rm *.gz*