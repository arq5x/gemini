####################################################################
# 1. Test accessing individual genotype alleles
####################################################################
echo "    genotypes.t01...\c"
echo "./.	./.	./.	./.
G/G	G/G	G/G	G/G
CCT/CCT	CCT/CCT	CCT/C	CCT/CCT
T/C	T/C	T/T	T/T
./.	./.	./.	./.
./.	./.	G/G	G/G
T/T	T/T	T/T	T/T
./.	./.	A/G	A/G
A/A	A/T	A/A	A/A
./.	G/G	G/G	G/G" > exp
gemini query -q "select gts.1094PC0005, gts.1094PC0009, \
				gts.1094PC0012, gts.1094PC0013 \
				from variants" test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 2. Test accessing individual genotype types
####################################################################
echo "    genotypes.t02...\c"
echo "2	2	2	2
0	0	0	0
0	0	1	0
1	1	0	0
2	2	2	2
2	2	3	3
0	0	0	0
2	2	1	1
0	1	0	0
2	0	0	0" > exp
gemini query -q "select gt_types.1094PC0005, gt_types.1094PC0009, \
	                    gt_types.1094PC0012, gt_types.1094PC0013 \
	             from variants" test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 3. Test filtering individual genotype types
####################################################################
echo "    genotypes.t03...\c"
echo "0	0	1	0
2	2	1	1" > exp
gemini query -q "select gt_types.1094PC0005, gt_types.1094PC0009, \
	                    gt_types.1094PC0012, gt_types.1094PC0013 \
	             from variants" \
			 --gt-filter "gt_types.1094PC0012 == HET" \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 4. Test a more complex genotype filter
####################################################################
echo "    genotypes.t04...\c"
echo "chrom	end	ref	alt	gt_types.1094PC0005	gt_types.1094PC0009	gt_types.1094PC0012	gt_types.1094PC0013
chr1	30869	CCT	C	0	0	1	0
chr1	30895	T	C	1	1	0	0
chr1	69511	A	G	2	2	1	1" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_types.1094PC0005, gt_types.1094PC0009, \
	                    gt_types.1094PC0012, gt_types.1094PC0013 \
	             from variants" \
			 --gt-filter "(gt_types.1094PC0012 == HET or \
						   gt_types.1094PC0005 == HET)" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 5.  Test accessing individual genotype depths
####################################################################
echo "    genotypes.t05...\c"
echo "chrom	end	ref	alt	gt_depths.1094PC0005	gt_depths.1094PC0009	gt_depths.1094PC0012	gt_depths.1094PC0013
chr1	30548	T	G	-1	-1	-1	-1
chr1	30860	G	C	7	2	6	4
chr1	30869	CCT	C	8	3	6	5
chr1	30895	T	C	8	3	6	5
chr1	30923	G	T	-1	-1	-1	-1
chr1	69270	A	G	-1	-1	3	6
chr1	69428	T	G	2	79	87	107
chr1	69511	A	G	-1	-1	6	4
chr1	69761	A	T	1	7	12	9
chr1	69871	G	A	-1	4	2	2" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_depths.1094PC0005, gt_depths.1094PC0009, \
	                    gt_depths.1094PC0012, gt_depths.1094PC0013 \
	             from variants" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

