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

# make sure a missing field works (1 recoord has not BaseQRankSum
gemini annotate -f test.anno.vcf.gz -o uniq_list -c BaseQRankSum -t float test.snpeff.vcf.db
gemini query -q "select chrom, end, BaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chr1	30548	0.0
chr1	30860	2.815
chr1	30869	0.9
chr1	30895	-0.23
chr1	30923	None
chr1	69270	1.675
chr1	69428	-0.638
chr1	69511	-1.592
chr1	69761	1.333
chr1	69871	-0.646" > exp

check obs exp test-annotate-vcf.t02



gemini annotate -f test.anno.vcf.gz -o uniq_list -e BaseQRankSum -c MyBaseQRankSum -t float test.snpeff.vcf.db
gemini query --header -q "select chrom, end, MyBaseQRankSum, BaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chrom	end	MyBaseQRankSum	BaseQRankSum
chr1	30548	0.0	0.0
chr1	30860	2.815	2.815
chr1	30869	0.9	0.9
chr1	30895	-0.23	-0.23
chr1	30923	None	None
chr1	69270	1.675	1.675
chr1	69428	-0.638	-0.638
chr1	69511	-1.592	-1.592
chr1	69761	1.333	1.333
chr1	69871	-0.646	-0.646" > exp

check obs exp test-annotate-vcf.t03

# check first

gemini annotate -f test.anno.vcf.gz -o first -e BaseQRankSum -c MyBaseQRankSum -t float test.snpeff.vcf.db
gemini query --header -q "select chrom, end, MyBaseQRankSum, BaseQRankSum from variants" test.snpeff.vcf.db > obs
echo "chrom	end	MyBaseQRankSum	BaseQRankSum
chr1	30548	0.0	0.0
chr1	30860	2.815	2.815
chr1	30869	0.9	0.9
chr1	30895	-0.23	-0.23
chr1	30923	None	None
chr1	69270	1.675	1.675
chr1	69428	-0.638	-0.638
chr1	69511	-1.592	-1.592
chr1	69761	1.333	1.333
chr1	69871	-0.646	-0.646" > exp

check obs exp test-annotate-vcf.t04

####################################
# create a missing record. and see -a boolean
####################################
zgrep -vw 69270 test.anno.vcf.gz | bgzip -c > t_.vcf.gz
tabix t_.vcf.gz

gemini annotate -f t_.vcf.gz -a boolean test.snpeff.vcf.db -c BaseQRankSum

gemini query --header -q "select chrom, end, MyBaseQRankSum, BaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chrom	end	MyBaseQRankSum	BaseQRankSum
chr1	30548	0.0	1.0
chr1	30860	2.815	1.0
chr1	30869	0.9	1.0
chr1	30895	-0.23	1.0
chr1	30923	None	1.0
chr1	69270	1.675	0.0
chr1	69428	-0.638	1.0
chr1	69511	-1.592	1.0
chr1	69761	1.333	1.0
chr1	69871	-0.646	1.0" > exp

check obs exp test-annotate-vcf.t05

####################################
# create a missing record. and see -a count
####################################

gemini annotate -f t_.vcf.gz -a count test.snpeff.vcf.db -c CountBaseQRankSum
gemini query --header -q "select chrom, end, MyBaseQRankSum, BaseQRankSum, CountBaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chrom	end	MyBaseQRankSum	BaseQRankSum	CountBaseQRankSum
chr1	30548	0.0	1.0	1
chr1	30860	2.815	1.0	1
chr1	30869	0.9	1.0	1
chr1	30895	-0.23	1.0	1
chr1	30923	None	1.0	1
chr1	69270	1.675	0.0	0
chr1	69428	-0.638	1.0	1
chr1	69511	-1.592	1.0	1
chr1	69761	1.333	1.0	1
chr1	69871	-0.646	1.0	1" > exp

check obs exp test-annotate-vcf.t06

############################################
# TODO: check mean, etc. with missing values.
############################################

gunzip t_.vcf.gz
grep -w 69871 t_.vcf > t_
cat t_ >> t_.vcf; rm t_
bgzip -f t_.vcf
tabix t_.vcf.gz

gemini annotate -f t_.vcf.gz -o sum -e BaseQRankSum -c SumBaseQRankSum -t float test.snpeff.vcf.db
gemini query --header -q "select chrom, end, BaseQRankSum, SumBaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chrom	end	BaseQRankSum	SumBaseQRankSum
chr1	30548	1.0	None
chr1	30860	1.0	2.815
chr1	30869	1.0	0.9
chr1	30895	1.0	-0.23
chr1	30923	1.0	None
chr1	69270	0.0	None
chr1	69428	1.0	-0.638
chr1	69511	1.0	-1.592
chr1	69761	1.0	1.333
chr1	69871	1.0	-1.292" > exp # note how the last number here is double of -0.646 becaus ewe included that last record twice.

check obs exp test-annotate-vcf.t06


gemini annotate -f t_.vcf.gz -a count test.snpeff.vcf.db -c CountBaseQRankSum
gemini query --header -q "select chrom, end, CountBaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chrom	end	CountBaseQRankSum
chr1	30548	1
chr1	30860	1
chr1	30869	1
chr1	30895	1
chr1	30923	1
chr1	69270	0
chr1	69428	1
chr1	69511	1
chr1	69761	1
chr1	69871	2" > exp # note how the lat one has value of 2.

check obs exp test-annotate-vcf.t07

############################
# Check Missing Values
############################

# add a 3rd value to the last row
zless t_.vcf.gz > t2_.vcf
grep -w 69871 t2_.vcf | perl -pe 's/BaseQRankSum=([^;]+);//' > t_
cat t_ >> t2_.vcf; rm t_
bgzip -f t2_.vcf
tabix t2_.vcf.gz



gemini annotate -f t2_.vcf.gz -t float -e BaseQRankSum -o mean test.snpeff.vcf.db -c Missing
gemini query --header -q "select chrom, end, Missing from variants" test.snpeff.vcf.db > obs

echo "chrom	end	Missing
chr1	30548	0
chr1	30860	2.815
chr1	30869	0.9
chr1	30895	-0.23
chr1	30923	None
chr1	69270	None
chr1	69428	-0.638
chr1	69511	-1.592
chr1	69761	1.333
chr1	69871	-0.646" > exp

check obs exp test-annotate-vcf.t08

### NOTE!!!
# that sum is 0 since np.sum([]) == 0 and mean is nan since np.mean([]) == nan
###

gemini annotate -f t2_.vcf.gz -t float -e BaseQRankSum -o sum test.snpeff.vcf.db -c Missing
gemini query --header -q "select chrom, end, Missing from variants" test.snpeff.vcf.db > obs

echo "chrom	end	Missing
chr1	30548	0
chr1	30860	2.815
chr1	30869	0.9
chr1	30895	-0.23
chr1	30923	0
chr1	69270	None
chr1	69428	-0.638
chr1	69511	-1.592
chr1	69761	1.333
chr1	69871	-1.292" > exp

check obs exp test-annotate-vcf.t09

## test multiple columns
gemini annotate -f test.anno.vcf.gz -t float,float -e AN,BaseQRankSum -o mean,list test.snpeff.vcf.db 
gemini query --header -q "select chrom, end, AN,BaseQRankSum from variants" test.snpeff.vcf.db > obs

echo "chrom	end	AN	BaseQRankSum
chr1	30548	14.0	0.0
chr1	30860	52.0	2.815
chr1	30869	56.0	0.9
chr1	30895	58.0	-0.23
chr1	30923	30.0	None
chr1	69270	56.0	1.675
chr1	69428	106.0	-0.638
chr1	69511	90.0	-1.592
chr1	69761	88.0	1.333
chr1	69871	80.0	-0.646" > exp

check obs exp test-annotate-vcf.t10


## test multiple columns with alternate names
gemini annotate -f test.anno.vcf.gz -t float,float -e AN,BaseQRankSum -c abc,def -o mean,list test.snpeff.vcf.db 
gemini query --header -q "select chrom, end, abc, def from variants" test.snpeff.vcf.db > obs

echo "chrom	end	abc	def
chr1	30548	14.0	0.0
chr1	30860	52.0	2.815
chr1	30869	56.0	0.9
chr1	30895	58.0	-0.23
chr1	30923	30.0	None
chr1	69270	56.0	1.675
chr1	69428	106.0	-0.638
chr1	69511	90.0	-1.592
chr1	69761	88.0	1.333
chr1	69871	80.0	-0.646" > exp

check obs exp test-annotate-vcf.t11
