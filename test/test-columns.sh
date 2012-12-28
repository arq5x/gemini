####################################################################
# 1. Make sure we are representing the correct REF and ALT alleles.
####################################################################
echo "    columns.t01...\c"
echo "chr1	30547	30548	T	G
chr1	30859	30860	G	C
chr1	30866	30869	CCT	C
chr1	30894	30895	T	C
chr1	30922	30923	G	T
chr1	69269	69270	A	G
chr1	69427	69428	T	G
chr1	69510	69511	A	G
chr1	69760	69761	A	T
chr1	69870	69871	G	A" > exp
gemini query -q "select chrom, start, end, ref, alt from variants" \
	   test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 2. Test variant type and sub_type
####################################################################
echo "    columns.t02...\c"
echo "snp	tv
snp	tv
indel	del
snp	ts
snp	tv
snp	ts
snp	tv
snp	ts
snp	tv
snp	ts" > exp
gemini query -q "select type, sub_type from variants" \
	   test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 3. Test depth
####################################################################
echo "    columns.t03...\c"
echo "7
61
65
63
17
100
5013
3853
574
154" > exp
gemini query -q "select depth from variants" \
	   test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 4. Test qual
####################################################################
echo "    columns.t04...\c"
echo "50.09
54.3
49.48
51.79
601.49
2971.02
938.56
101944.44
1026.71
82.97" > exp
gemini query -q "select qual from variants" \
	   test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 5. Test filter
####################################################################
echo "    columns.t05...\c"
echo "PASS
PASS
PASS
PASS
PASS
PASS
PASS
PASS
PASS
PASS" > exp
gemini query -q "select filter from variants" \
	   test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 6. Test genotype counts
####################################################################
echo "    columns.t06...\c"
echo "0	0	0	4
4	0	0	0
3	1	0	0
2	2	0	0
0	0	0	4
0	0	2	2
4	0	0	0
0	2	0	2
3	1	0	0
3	0	0	1" > exp
gemini query -q "select num_hom_ref, num_het, num_hom_alt, num_unknown \
	             from variants" test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 7. Test call rate
####################################################################
echo "    columns.t07...\c"
echo "0.0
1.0
1.0
1.0
0.0
0.5
1.0
0.5
1.0
0.75" > exp
gemini query -q "select call_rate \
	             from variants" test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 8. Test alternate allele frequency 
####################################################################
echo "    columns.t08...\c"
echo "0.0
0.0
0.125
0.25
0.0
1.0
0.0
0.5
0.125
0.0" > exp
gemini query -q "select aaf \
	             from variants" test.snpEff.vcf.db \
       > obs
check obs exp
rm obs exp





