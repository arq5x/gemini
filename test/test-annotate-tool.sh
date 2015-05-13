###########################################################################################
#1. Test annotating variants using the "boolean" function
###########################################################################################
source ./check.sh

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548
chr1	30922	30923" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -c anno -a boolean test.snpeff.vcf.db

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
check obs exp annotate-tool.t1
rm obs exp
rm *.gz*


###########################################################################################
#2. Test annotating variants using the "count" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548
chr1	30920	30925
chr1	30922	30923" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -c anno2 -a count test.snpeff.vcf.db

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
check obs exp annotate-tool.t2
rm obs exp
rm *.gz*


###########################################################################################
#3. Test annotating variants using the "extract" function
#   while extacting just one column
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a
chr1	30920	30925	b
chr1	30922	30923	c" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -c anno3 -a extract -e 4 -t text -o list  test.snpeff.vcf.db

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
check obs exp annotate-tool.t3
rm obs exp
rm *.gz*


###########################################################################################
#4. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno4,anno5 -e 4,5 -t text,float -o list,mean  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	b,c	0.3
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno4, anno5 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t4
rm obs exp
rm *.gz*


###########################################################################################
#5. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	c	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno6,anno7 -e 4,5 -t text,float -o list,mean  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	b,c,c	0.4
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno6, anno7 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t5
rm obs exp
rm *.gz*


###########################################################################################
#6. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	c	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno8,anno9 -e 4,5 -t text,float -o uniq_list,mean  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	c,b	0.4
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno8, anno9 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t6
rm obs exp
rm *.gz*


###########################################################################################
#7. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	d	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno10,anno11 -e 4,5 -t text,float -o first,mean  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	b	0.4
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno10, anno11 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t7
rm obs exp
rm *.gz*

###########################################################################################
#8. Test annotating variants using the "extract" function
###########################################################################################
echo "    annotate-tool.t8...\c"

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	d	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno12,anno13 -e 4,5 -t text,float -o first,mean  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	b	0.4
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno12, anno13 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t8
rm obs exp
rm *.gz*

###########################################################################################
#9. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	d	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno14,anno15 -e 4,5 -t text,float -o first,last  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	b	0.6
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno14, anno15 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t9
rm obs exp
rm *.gz*


###########################################################################################
#10. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	d	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno16,anno17 -e 4,5 -t text,float -o last,first  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	d	0.2
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno16, anno17 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annoate-tool.t10
rm obs exp
rm *.gz*


###########################################################################################
#11. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	d	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno18,anno19 -e 4,5 -t text,float -o last,max  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	d	0.6
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno18, anno19 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t11
rm obs exp
rm *.gz*


###########################################################################################
#12. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.4
chr1	30922	30923	d	0.6" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno20,anno21 -e 4,5 -t text,float -o last,min  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	d	0.2
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno20, anno21 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t12
rm obs exp
rm *.gz*

###########################################################################################
#13. Test annotating variants using the "extract" function
###########################################################################################
echo "    annotate-tool.t13...\c"

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.3
chr1	30922	30923	d	0.4
chr1	30922	30923	e	0.5" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno22,anno23 -e 4,5 -t text,float -o last,median  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	e	0.35
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno22, anno23 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t13
rm obs exp
rm *.gz*


###########################################################################################
#14. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.3
chr1	30922	30923	d	0.3
chr1	30922	30923	e	0.5" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
gemini annotate -f anno.bed.gz -a extract -c anno24,anno25 -e 4,5 -t text,float -o last,mode  test.snpeff.vcf.db

echo "chr1	30548	a	0.23
chr1	30860	None	None
chr1	30869	None	None
chr1	30895	None	None
chr1	30923	e	0.3
chr1	69270	None	None
chr1	69428	None	None
chr1	69511	None	None
chr1	69761	None	None
chr1	69871	None	None" > exp

gemini query -q "select chrom, end, anno24, anno25 from variants" \
	test.snpeff.vcf.db > obs
check obs exp annotate-tool.t14
rm obs exp
rm *.gz*



###########################################################################################
#15. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.3
chr1	30922	30923	d	0.3
chr1	30922	30923	e	0.5" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
echo $'EXITING: The number of column names, numbers, types, and operations must match: [anno23], [4,5], [text,float], [last,mode]\n' > exp

gemini annotate -f anno.bed.gz -a extract -c anno23 -e 4,5 -t text,float -o last,mode  test.snpeff.vcf.db 2> obs

check obs exp annotate-tool.t15
rm obs exp
rm *.gz*


##########################################################################################
#16. Test annotating variants using the "extract" function
##########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.3
chr1	30922	30923	d	0.3
chr1	30922	30923	e	0.5" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
echo "usage: gemini annotate [-h] [-f ANNO_FILE] [-c COL_NAMES]
                       [-a {boolean,count,extract}] [-e COL_EXTRACTS]
                       [-t COL_TYPES] [-o COL_OPERATIONS] [--region-only]
                       db
gemini annotate: error: argument -a: invalid choice: 'distract' (choose from 'boolean', 'count', 'extract')" > exp

gemini annotate -f anno.bed.gz -a distract -c anno23,anno24 -e 4,5 -t text,float -o last,mode  test.snpeff.vcf.db 2> obs


check obs exp annotate-tool.t16
rm obs exp
rm *.gz*


###########################################################################################
#17. Test annotating variants using the "extract" function
###########################################################################################

# make a dunnmy TABIX'ed annotation file
echo "chr1	30547	30548	a	0.23
chr1	30920	30925	b	0.2
chr1	30922	30923	c	0.3
chr1	30922	30923	d	0.3
chr1	30922	30923	e	0.5" > anno.bed
bgzip anno.bed
tabix -p bed anno.bed.gz

# create a new column in the database using the new annotation
echo $'EXITING: Column operation [model] not supported.\n' > exp

gemini annotate -f anno.bed.gz -a extract -c anno27,anno28 -e 4,5 -t text,float -o last,model  test.snpeff.vcf.db 2> obs

check obs exp annotate-tool.t17
#rm obs exp
rm *.gz*
