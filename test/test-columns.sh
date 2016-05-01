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
	   test.snpeff.vcf.db \
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
	   test.snpeff.vcf.db \
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
	   test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 4. Test qual
####################################################################
echo "    columns.t04...\c"
echo "50.0900001526
54.2999992371
49.4799995422
51.7900009155
601.489990234
2971.02001953
938.559997559
101944.4375
1026.70996094
82.9700012207" > exp
gemini query -q "select qual from variants" \
	   test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 5. Test filter
####################################################################
echo "    columns.t05...\c"
echo "ABFilter
LowQual;QUALFilter
None
None
None
None
ABFilter;LowQual;QDFilter;QUALFilter;SBFilter
ABFilter;LowQual;QUALFilter
None" > exp
gemini query -q "select filter from variants" \
	   test2.snpeff.db \
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
	             from variants" test.snpeff.vcf.db \
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
	             from variants" test.snpeff.vcf.db \
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
	             from variants" test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

#####################################################################
# 9. Test allele_count
#####################################################################
echo "    columns.t09...\c"
echo "6
3
6
4
30
53
3
79
8
2" > exp
gemini query -q "select allele_count from variants" test.snpeff.vcf.db \
    > obs
    
check obs exp
rm obs exp

########################################################################
# 10. Test num_alleles
########################################################################
echo "    columns.t10...\c"
echo "14
52
56
58
30
56
106
90
88
80" > exp
gemini query -q "select num_alleles from variants" test.snpeff.vcf.db \
    > obs
check obs exp
rm obs exp

#########################################################################
# 11. Test qual_depth
#########################################################################
echo "    columns.t11...\c"
echo "16.7000007629
13.5699996948
3.80999994278
3.98000001907
35.3800010681
31.6100006104
18.3999996185
26.5100002289
21.3899993896
20.7399997711" > exp
gemini query -q "select qual_depth from variants" test.snpeff.vcf.db \
    > obs
check obs exp
rm obs exp

##########################################################################
# 12. Test rms_map_qual
###########################################################################
echo "    columns.t12...\c"
echo "29.0
36.25
36.0900001526
36.0900001526
35.6100006104
31.0599994659
29.7900009155
33.2400016785
32.1800003052
31.1000003815" > exp
gemini query -q "select rms_map_qual from variants" test.snpeff.vcf.db \
    > obs
check obs exp
rm obs exp

##############################################################################
#13. Test num map qual zero
##############################################################################
echo "    columns.t13...\c"
echo "0
0
0
0
0
0
0
0
0
0" > exp
gemini query -q "select num_mapq_zero from variants" test.snpeff.vcf.db \
 > obs
check obs exp
rm obs exp


##############################################################################
#14. Test the ID column
##############################################################################
echo "    columns.t14...\c"
echo "None
rs1234
None
rs567,rs89
None
foo
777" > exp
gemini query -q "select vcf_id from variants" test.vcf_id.snpeff.vcf.db \
 > obs
check obs exp
rm obs exp
