###################################################################
# 1. Test variants by sample
###################################################################
echo "    stat.t01...\c"
echo "sample	total
M10475	4
M10478	6
M10500	6
M128215	6" > exp
gemini stats --vars-by-sample test5.snpeff.db \
       > obs
check obs exp
rm obs exp

###################################################################
# 2. Test genotypes by sample
###################################################################
echo "    stat.t02...\c"
echo "sample	num_hom_ref	num_het	num_hom_alt	num_unknown	total
M10475	4	1	3	1	9
M10478	2	2	4	1	9
M10500	2	1	5	1	9
M128215	2	3	3	1	9" > exp
gemini stats --gts-by-sample test5.snpeff.db \
       > obs
check obs exp
rm obs exp

###################################################################
# 3. Test site freq. spectrum
###################################################################
echo "    stat.t03...\c"
echo "aaf	count
0.125	2
0.375	1
0.5	2
1.0	4" > exp
gemini stats --sfs test5.snpeff.db \
       > obs
check obs exp
rm obs exp

###################################################################
# 4. Test snp counts
###################################################################
echo "    stat.t04...\c"
echo "type	count
A->G	2
C->T	1
G->A	1
G->C	1
G->T	1" > exp
gemini stats --snp-counts test3.snpeff.db \
        > obs

check obs exp
rm obs exp

###################################################################
# 5. Test transition / transversion ratios
###################################################################
echo "    stat.t05...\c"
echo "ts	tv	ts/tv
4	5	0.8" > exp
gemini stats --tstv test.snpeff.vcf.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test transition / transversion ratios in coding region
###################################################################
echo "    stat.t06...\c"
echo "ts	tv	ts/tv
3	2	1.5" > exp
gemini stats --tstv-coding test.snpeff.vcf.db > obs
check obs exp
rm obs exp

###################################################################
# 7. Test transition / transversion ratios in noncoding region
###################################################################
echo "    stat.t07...\c"
echo "ts	tv	ts/tv
1	3	0.3333" > exp
gemini stats --tstv-noncoding test.snpeff.vcf.db > obs
check obs exp
rm obs exp


###################################################################
# 8. Test multi-dimensional scaling (mds)
###################################################################
echo "    stat.t08...\c"
echo "sample1	sample2	distance
M10500	M10500	0.0
M10475	M10475	0.0
M128215	M10500	2.5
M10478	M10478	0.0
M10475	M10500	2.0
M10500	M128215	2.5
M10475	M10478	1.25
M10500	M10475	2.0
M128215	M10478	1.7143
M10478	M10500	0.5714
M10475	M128215	0.5714
M10478	M128215	1.7143
M128215	M10475	0.5714
M128215	M128215	0.0
M10478	M10475	1.25
M10500	M10478	0.5714" > exp
gemini stats --mds test5.snpeff.db \
      > obs
check obs exp
rm obs exp

###################################################################
# 9. Test summarize
###################################################################
echo "    stat.t09...\c"
echo "sample	total	num_het	num_hom_alt	num_hom_ref
M10475	4	1	3	4
M128215	6	3	3	2
M10478	6	2	4	2
M10500	6	1	5	2" > exp
gemini stats --summarize "select * from variants" test5.snpeff.db > obs
check obs exp
rm obs exp
