###########################################################################################
#Test annotating variants with a VCF.
###########################################################################################
source ./check.sh

# make a dunnmy TABIX'ed annotation file
#BaseQRankSum
bgzip -c test.snpeff.vcf > test.anno.vcf.gz
tabix test.anno.vcf.gz

# create a new column in the database using the new annotation
gemini annotate -f test.anno.vcf.gz -o uniq_list -c AN -t float test.snpeff.vcf.db

gemini query -q "select chrom, end, AN from variants" test.snpeff.vcf.db > obs
echo "chr1	30548	14.0
chr1	30860	52.0
chr1	30869	56.0
chr1	30895	58.0
chr1	30923	30.0
chr1	69270	56.0
chr1	69428	106.0
chr1	69511	90.0
chr1	69761	88.0
chr1	69871	80.0" > exp

check obs exp test-annotate-vcf.t01
