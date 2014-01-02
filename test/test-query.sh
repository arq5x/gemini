check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

####################################################################
# 1. Test the samples table
####################################################################
echo "    query.t01...\c"
echo "1	None	1094PC0005	None	None	None	None
2	None	1094PC0009	None	None	None	None
3	None	1094PC0012	None	None	None	None
4	None	1094PC0013	None	None	None	None
5	None	1094PC0016	None	None	None	None
6	None	1094PC0017	None	None	None	None
7	None	1094PC0018	None	None	None	None
8	None	1094PC0019	None	None	None	None
9	None	1094PC0020	None	None	None	None
10	None	1094PC0021	None	None	None	None" > exp
gemini query -q "select * from samples limit 10" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 2. Test a basic query of the variants table
####################################################################
echo "    query.t02...\c"
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
gemini query -q "select chrom, start, end, ref, alt from variants limit 10" test.query.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 3. Test a basic query of the variants table with a where clause
####################################################################
echo "    query.t03...\c"
echo "chr1	1219381	1219382	C	G	SCNN1D
chr1	1219476	1219477	T	G	SCNN1D
chr1	1219486	1219487	T	G	SCNN1D
chr1	1219488	1219489	A	G	SCNN1D
chr1	1219494	1219496	GT	G	SCNN1D
chr1	1219502	1219505	GTT	G	SCNN1D
chr1	1219507	1219511	GTGA	G	SCNN1D
chr1	1219521	1219524	GTC	G	SCNN1D
chr1	1219533	1219536	GTT	G	SCNN1D
chr1	1219555	1219558	GTT	G	SCNN1D" > exp
gemini query -q "select chrom, start, end, ref, alt, gene \
                 from variants \
                 where gene == 'SCNN1D' limit 10" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 4. Test a query of the variants table with a where clause
#    and a request for a sample's genotype
####################################################################
echo "    query.t04...\c"
echo "chr1	1219381	1219382	C	G	SCNN1D	C/C
chr1	1219476	1219477	T	G	SCNN1D	T/T
chr1	1219486	1219487	T	G	SCNN1D	T/T
chr1	1219488	1219489	A	G	SCNN1D	A/A
chr1	1219494	1219496	GT	G	SCNN1D	GT/GT
chr1	1219502	1219505	GTT	G	SCNN1D	./.
chr1	1219507	1219511	GTGA	G	SCNN1D	./.
chr1	1219521	1219524	GTC	G	SCNN1D	./.
chr1	1219533	1219536	GTT	G	SCNN1D	./.
chr1	1219555	1219558	GTT	G	SCNN1D	./." > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gts.1094PC0018 \
                 from variants \
                 where gene == 'SCNN1D' limit 10" test.query.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 5. Test a query of the variants table with a where clause
#    and a request for a sample's genotype and type
####################################################################
echo "    query.t05...\c"
echo "chr1	1219381	1219382	C	G	SCNN1D	C/C	0
chr1	1219476	1219477	T	G	SCNN1D	T/T	0
chr1	1219486	1219487	T	G	SCNN1D	T/T	0
chr1	1219488	1219489	A	G	SCNN1D	A/A	0
chr1	1219494	1219496	GT	G	SCNN1D	GT/GT	0" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gts.1094PC0018, gt_types.1094PC0018 \
                 from variants \
                 where gene == 'SCNN1D' limit 5" test.query.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 6. Test a query of the variants table with a
# request for all sample genotypes
####################################################################
echo "    query.t06...\c"
echo "chr1	30547	30548	T	G	FAM138A	./.,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,G/G,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,T/T,./.,G/G,./.,G/G,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,./.
chr1	30859	30860	G	C	FAM138A	G/G,G/G,G/G,G/G,./.,./.,G/G,./.,G/G,G/G,G/G,./.,./.,G/G,G/G,./.,./.,./.,./.,./.,G/G,./.,./.,G/G,./.,G/G,G/G,G/G,G/G,C/C,G/G,./.,G/G,./.,./.,./.,./.,./.,./.,./.,./.,./.,G/C,./.,G/G,./.,G/G,./.,./.,G/G,G/G,G/G,G/G,./.,./.,./.,./.,./.,./.,./.
chr1	30866	30869	CCT	C	FAM138A	CCT/CCT,CCT/CCT,CCT/C,CCT/CCT,./.,./.,CCT/CCT,./.,CCT/CCT,CCT/C,CCT/CCT,./.,./.,CCT/CCT,CCT/CCT,./.,./.,./.,./.,./.,CCT/CCT,./.,./.,CCT/C,./.,CCT/CCT,C/C,CCT/CCT,CCT/CCT,CCT/CCT,CCT/CCT,./.,CCT/CCT,./.,./.,./.,./.,CCT/CCT,./.,./.,./.,./.,CCT/C,./.,CCT/CCT,./.,CCT/CCT,CCT/CCT,./.,CCT/CCT,CCT/CCT,CCT/CCT,CCT/CCT,./.,./.,./.,./.,./.,./.,./.
chr1	30894	30895	T	C	FAM138A	T/C,T/C,T/T,T/T,./.,./.,T/T,./.,T/T,T/T,T/T,./.,T/T,T/T,T/T,./.,./.,./.,T/T,./.,T/T,./.,./.,T/T,./.,./.,T/T,T/T,T/T,C/C,T/T,./.,T/T,./.,./.,./.,./.,T/T,./.,./.,./.,./.,T/T,./.,T/T,./.,T/T,T/T,./.,./.,T/T,T/T,T/T,T/T,./.,./.,./.,./.,./.,./.
chr1	30922	30923	G	T	FAM138A	./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,T/T,./.,T/T,./.,./.,./.,./.,T/T,T/T,T/T,T/T,./.,T/T,./.,T/T,./.,./.,./.,./.,T/T,T/T,./.,./.,./.,./.,./.,T/T,./.,T/T,T/T,./.,./.,./.,./.,T/T,T/T,./.,./.,./.,./.,./.,./." > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gts \
                 from variants \
                 limit 5" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 7. Test a query of the variants table with a
# request for all sample genotype types
####################################################################
echo "    query.t07...\c"
echo "chr1	30547	30548	T	G	FAM138A	2,2,2,2,2,2,0,2,2,2,2,2,2,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,0,2,3,2,3,2,2,2,2,2,0,2,2,2,2,2,2,2
chr1	30859	30860	G	C	FAM138A	0,0,0,0,2,2,0,2,0,0,0,2,2,0,0,2,2,2,2,2,0,2,2,0,2,0,0,0,0,3,0,2,0,2,2,2,2,2,2,2,2,2,1,2,0,2,0,2,2,0,0,0,0,2,2,2,2,2,2,2
chr1	30866	30869	CCT	C	FAM138A	0,0,1,0,2,2,0,2,0,1,0,2,2,0,0,2,2,2,2,2,0,2,2,1,2,0,3,0,0,0,0,2,0,2,2,2,2,0,2,2,2,2,1,2,0,2,0,0,2,0,0,0,0,2,2,2,2,2,2,2
chr1	30894	30895	T	C	FAM138A	1,1,0,0,2,2,0,2,0,0,0,2,0,0,0,2,2,2,0,2,0,2,2,0,2,2,0,0,0,3,0,2,0,2,2,2,2,0,2,2,2,2,0,2,0,2,0,0,2,2,0,0,0,0,2,2,2,2,2,2
chr1	30922	30923	G	T	FAM138A	2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,3,2,2,2,2,3,3,3,3,2,3,2,3,2,2,2,2,3,3,2,2,2,2,2,3,2,3,3,2,2,2,2,3,3,2,2,2,2,2,2" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gt_types \
                 from variants \
                 limit 5" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 8. Test a query of the variants table with a
# request for all sample genotype phases
####################################################################
echo "    query.t08...\c"
echo "chr1	30547	30548	T	G	FAM138A	False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False
chr1	30859	30860	G	C	FAM138A	False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False
chr1	30866	30869	CCT	C	FAM138A	False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False
chr1	30894	30895	T	C	FAM138A	False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False
chr1	30922	30923	G	T	FAM138A	False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False,False" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gt_phases \
                 from variants \
                 limit 5" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 9. Test a query of the variants table with a
# request for all sample genotype phases
####################################################################
echo "    query.t09...\c"
echo "chr1	30547	30548	T	G	FAM138A	-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1
chr1	30859	30860	G	C	FAM138A	7,2,6,4,-1,-1,1,-1,3,2,1,-1,-1,1,2,-1,-1,-1,-1,-1,1,-1,-1,2,-1,1,1,3,2,2,3,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,-1,1,-1,2,-1,-1,1,1,1,2,-1,-1,-1,-1,-1,-1,-1
chr1	30866	30869	CCT	C	FAM138A	8,3,6,5,-1,-1,1,-1,3,2,1,-1,-1,1,2,-1,-1,-1,-1,-1,1,-1,-1,2,-1,1,1,3,2,2,3,-1,2,-1,-1,-1,-1,1,-1,-1,-1,-1,2,-1,1,-1,2,1,-1,1,1,1,2,-1,-1,-1,-1,-1,-1,-1
chr1	30894	30895	T	C	FAM138A	8,3,6,5,-1,-1,1,-1,3,2,1,-1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,2,-1,-1,1,3,2,2,2,-1,2,-1,-1,-1,-1,1,-1,-1,-1,-1,2,-1,1,-1,2,1,-1,-1,1,1,2,1,-1,-1,-1,-1,-1,-1
chr1	30922	30923	G	T	FAM138A	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1,-1,-1,1,1,1,1,-1,1,-1,2,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,2,1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gt_depths \
                 from variants \
                 limit 5" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 10. Test a query of the variants table with a select *
####################################################################
echo "    query.t10...\c"
echo "chr1	30547	30548	None	1	1	T	G	50.09	None	snp	tv	0.116666666667	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	4	0	3	53	0.428571428571	0.00815097309095	None	0.527472527473	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	7	None	29.0	0	0	14	0.0	0.0	16.7	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown
chr1	30859	30860	None	2	1	G	C	54.3	None	snp	tv	0.433333333333	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	24	1	1	34	0.0576923076923	0.000983220876781	0.646258503401	0.110859728507	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	61	None	36.25	0	0	52	0.0	0.6396	13.57	3	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown
chr1	30866	30869	None	3	1	CCT	C	49.48	None	indel	del	0.466666666667	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	23	4	1	32	0.107142857143	0.180078296859	0.253333333333	0.194805194805	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	65	None	36.09	0	0	56	None	12.2055	3.81	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown
chr1	30894	30895	None	4	1	T	C	51.79	None	snp	ts	0.483333333333	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	26	2	1	31	0.0689655172414	0.0126621812074	0.462962962963	0.130671506352	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	63	None	36.09	0	0	58	0.0	0.2091	3.98	4	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown
chr1	30922	30923	None	5	1	G	T	601.49	None	snp	tv	0.25	1	rs140337953	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	0	0	15	45	1.0	1	None	0	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	17	None	35.61	0	0	30	0.0	0.0	35.38	30	None	None	None	None	0	None	None	None	0	1	0.8	0.89	0.48	0.73	0.73	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown" > exp
gemini query -q "select * \
                 from variants \
                 limit 5" test.query.db \
       > obs
check obs exp
rm obs exp


##############################################################################
# 11. Test a query of the variants table with a select * and a genotype column
###############################################################################
echo "    query.t11...\c"
echo "chr1	30547	30548	None	1	1	T	G	50.09	None	snp	tv	0.116666666667	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	4	0	3	53	0.428571428571	0.00815097309095	None	0.527472527473	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	7	None	29.0	0	0	14	0.0	0.0	16.7	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	T/T
chr1	30859	30860	None	2	1	G	C	54.3	None	snp	tv	0.433333333333	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	24	1	1	34	0.0576923076923	0.000983220876781	0.646258503401	0.110859728507	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	61	None	36.25	0	0	52	0.0	0.6396	13.57	3	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	G/G
chr1	30866	30869	None	3	1	CCT	C	49.48	None	indel	del	0.466666666667	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	23	4	1	32	0.107142857143	0.180078296859	0.253333333333	0.194805194805	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	65	None	36.09	0	0	56	None	12.2055	3.81	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	CCT/CCT
chr1	30894	30895	None	4	1	T	C	51.79	None	snp	ts	0.483333333333	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	26	2	1	31	0.0689655172414	0.0126621812074	0.462962962963	0.130671506352	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	63	None	36.09	0	0	58	0.0	0.2091	3.98	4	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	T/T
chr1	30922	30923	None	5	1	G	T	601.49	None	snp	tv	0.25	1	rs140337953	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	0	0	15	45	1.0	1	None	0	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	17	None	35.61	0	0	30	0.0	0.0	35.38	30	None	None	None	None	0	None	None	None	0	1	0.8	0.89	0.48	0.73	0.73	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	./." > exp
gemini query -q "select *, gts.1094PC0018 \
                  from variants \
                  limit 5" test.query.db \
        > obs
check obs exp
rm obs exp


##############################################################################
# 12. Test a query of the variants table with a select * and the full genotype column
###############################################################################
echo "    query.t12...\c"
echo "chr1	30547	30548	None	1	1	T	G	50.09	None	snp	tv	0.116666666667	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.3	None	0	1	0	None	None	4	0	3	53	0.428571428571	0.00815097309095	None	0.527472527473	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	8	protein_coding	downstream	LOW	None	None	None	None	None	None	None	7	None	29.0	0	0	14	0.0	0.0	16.7	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	None	./.,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,G/G,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,T/T,./.,G/G,./.,G/G,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,./.
chr1	30859	30860	None	2	1	G	C	54.3	None	snp	tv	0.433333333333	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.3	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	24	1	1	34	0.0576923076923	0.000983220876781	0.646258503401	0.110859728507	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	61	None	36.25	0	0	52	0.0	0.6396	13.57	3	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	None	G/G,G/G,G/G,G/G,./.,./.,G/G,./.,G/G,G/G,G/G,./.,./.,G/G,G/G,./.,./.,./.,./.,./.,G/G,./.,./.,G/G,./.,G/G,G/G,G/G,G/G,C/C,G/G,./.,G/G,./.,./.,./.,./.,./.,./.,./.,./.,./.,G/C,./.,G/G,./.,G/G,./.,./.,G/G,G/G,G/G,G/G,./.,./.,./.,./.,./.,./.,./.
chr1	30866	30869	None	3	1	CCT	C	49.48	None	indel	del	0.466666666667	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.3	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	23	4	1	32	0.107142857143	0.180078296859	0.253333333333	0.194805194805	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	65	None	36.09	0	0	56	None	12.2055	3.81	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	None	CCT/CCT,CCT/CCT,CCT/C,CCT/CCT,./.,./.,CCT/CCT,./.,CCT/CCT,CCT/C,CCT/CCT,./.,./.,CCT/CCT,CCT/CCT,./.,./.,./.,./.,./.,CCT/CCT,./.,./.,CCT/C,./.,CCT/CCT,C/C,CCT/CCT,CCT/CCT,CCT/CCT,CCT/CCT,./.,CCT/CCT,./.,./.,./.,./.,CCT/CCT,./.,./.,./.,./.,CCT/C,./.,CCT/CCT,./.,CCT/CCT,CCT/CCT,./.,CCT/CCT,CCT/CCT,CCT/CCT,CCT/CCT,./.,./.,./.,./.,./.,./.,./.
chr1	30894	30895	None	4	1	T	C	51.79	None	snp	ts	0.483333333333	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.3	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	26	2	1	31	0.0689655172414	0.0126621812074	0.462962962963	0.130671506352	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	63	None	36.09	0	0	58	0.0	0.2091	3.98	4	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	None	T/C,T/C,T/T,T/T,./.,./.,T/T,./.,T/T,T/T,T/T,./.,T/T,T/T,T/T,./.,./.,./.,T/T,./.,T/T,./.,./.,T/T,./.,./.,T/T,T/T,T/T,C/C,T/T,./.,T/T,./.,./.,./.,./.,T/T,./.,./.,./.,./.,T/T,./.,T/T,./.,T/T,T/T,./.,./.,T/T,T/T,T/T,T/T,./.,./.,./.,./.,./.,./.
chr1	30922	30923	None	5	1	G	T	601.49	None	snp	tv	0.25	1	rs140337953	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.3	Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A	0	1	0	None	None	0	0	15	45	1.0	1	None	0	2.981822	FAM138A	ENST00000417324	0	0	0	None	None	None	85	protein_coding	downstream	LOW	None	None	None	None	None	None	None	17	None	35.61	0	0	30	0.0	0.0	35.38	30	None	None	None	None	0	None	None	None	0	1	0.8	0.89	0.48	0.73	0.73	None	None	None	None	0	None	None	None	R	R	R	R	R	unknown	None	./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,T/T,./.,T/T,./.,./.,./.,./.,T/T,T/T,T/T,T/T,./.,T/T,./.,T/T,./.,./.,./.,./.,T/T,T/T,./.,./.,./.,./.,./.,T/T,./.,T/T,T/T,./.,./.,./.,./.,T/T,T/T,./.,./.,./.,./.,./.,./." > exp
gemini query -q "select *, gts \
                  from variants \
                  limit 5" test.query.db > obs
check obs exp
rm obs exp

####################################################################
# 13. Test a query of the variants table with a where clause
#     and a genotype filter
####################################################################
echo "    query.t13...\c"
echo "chr1	1219381	1219382	C	G	SCNN1D	C/C	0
chr1	1219476	1219477	T	G	SCNN1D	T/T	0
chr1	1219486	1219487	T	G	SCNN1D	T/T	0
chr1	1219488	1219489	A	G	SCNN1D	A/A	0
chr1	1219494	1219496	GT	G	SCNN1D	GT/GT	0" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gts.1094PC0018, gt_types.1094PC0018 \
                 from variants \
                 where gene == 'SCNN1D' limit 5" \
             --gt-filter "gt_types.1094PC0018 != HET" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 14. Test a query of the variants table with a where clause
#     and a more complex genotype filter
####################################################################
echo "    query.t14...\c"
echo "chr1	1219381	1219382	C	G	SCNN1D	C/C	C/C
chr1	1219476	1219477	T	G	SCNN1D	T/T	T/T
chr1	1219486	1219487	T	G	SCNN1D	T/T	T/T
chr1	1219488	1219489	A	G	SCNN1D	A/A	A/A" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gts.1094PC0018, gts.1094PC0019 \
                 from variants \
                 where gene == 'SCNN1D' limit 5" \
             --gt-filter "gt_types.1094PC0018 == HET or gt_types.1094PC0019 == HOM_REF" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 15. Test a query of the variants table with a more complex genotype filter
####################################################################
echo "    query.t15...\c"
echo "chr1	865218	865219	G	A	SAMD11	G/A	G/G
chr1	874948	874949	G	GAC	SAMD11	G/GAC	G/G
chr1	880750	880751	A	G	SAMD11	A/G	A/A
chr1	886408	886409	G	C	NOC2L	G/C	G/G
chr1	891897	891898	T	C	NOC2L	T/C	T/T
chr1	892470	892471	G	A	NOC2L	G/A	G/G
chr1	898602	898603	C	G	KLHL17	C/G	C/C
chr1	906301	906302	C	T	C1orf170	C/T	C/C
chr1	908822	908823	G	A	C1orf170	G/A	G/G
chr1	909308	909309	T	C	PLEKHN1	T/C	T/T
chr1	909418	909419	C	T	C1orf170	C/T	C/C
chr1	934796	934797	T	G	HES4	T/G	T/T
chr1	970559	970563	GGGT	G	AGRN	GGGT/G	GGGT/GGGT
chr1	970561	970563	GT	G	AGRN	GT/G	GT/GT
chr1	978761	978762	G	A	AGRN	G/A	G/G
chr1	978856	978857	T	G	AGRN	T/G	T/T
chr1	979593	979594	C	T	AGRN	C/T	C/C
chr1	982843	982844	G	C	AGRN	G/C	G/G
chr1	985444	985445	G	GT	AGRN	G/GT	G/G
chr1	985445	985446	G	T	AGRN	G/T	G/G
chr1	985446	985447	G	T	AGRN	G/T	G/G
chr1	985661	985662	C	T	AGRN	C/T	C/C
chr1	986884	986885	T	G	AGRN	T/G	T/T
chr1	987310	987311	T	C	AGRN	T/C	T/T
chr1	1119542	1119543	G	C	TTLL10	G/C	G/G
chr1	1158325	1158326	G	A	SDF4	G/A	G/G
chr1	1158356	1158357	A	G	SDF4	A/G	A/A
chr1	1158440	1158443	GCA	G	SDF4	GCA/G	GCA/GCA
chr1	1158533	1158534	G	GAC	SDF4	G/GAC	G/G
chr1	1158561	1158564	AAC	A	SDF4	AAC/A	AAC/AAC
chr1	1158566	1158567	A	G	SDF4	A/G	A/A
chr1	1158947	1158948	A	G	SDF4	A/G	A/A
chr1	1158972	1158973	A	T	SDF4	A/T	A/A
chr1	1158974	1158975	A	C	SDF4	A/C	A/A
chr1	1159484	1159485	C	T	SDF4	C/T	C/C
chr1	1163803	1163804	C	T	SDF4	C/T	C/C
chr1	1179415	1179416	A	C	FAM132A	A/C	A/A
chr1	1181371	1181372	C	T	FAM132A	C/T	C/C
chr1	1192771	1192773	CA	C	UBE2J2	CA/C	CA/CA" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gts.1094PC0018, gts.1094PC0019 \
                 from variants" \
             --gt-filter "gt_types.1094PC0018 == HET and gt_types.1094PC0019 == HOM_REF" test.query.db \
       > obs
check obs exp
rm obs exp

#########################################################################
# 16. Test a basic query of the variants table with show-variant-samples
#########################################################################
echo "    query.t16...\c"
echo "chrom	start	end	ref	alt	variant_samples	HET_samples	HOM_ALT_samples
chr1	30547	30548	T	G	1478PC0016,1719PC0007,1719PC0009		1478PC0016,1719PC0007,1719PC0009
chr1	30859	30860	G	C	1478PC0017B,1719PC0005	1719PC0005	1478PC0017B
chr1	30866	30869	CCT	C	1094PC0012,1094PC0021,1478PC0011,1478PC0014B,1719PC0005	1094PC0012,1094PC0021,1478PC0011,1719PC0005	1478PC0014B
chr1	30894	30895	T	C	1094PC0005,1094PC0009,1478PC0017B	1094PC0005,1094PC0009	1478PC0017B
chr1	30922	30923	G	T	1478PC0006B,1478PC0008B,1478PC0013B,1478PC0014B,1478PC0015B,1478PC0016,1478PC0018,1478PC0020,1478PC0025,1719PC0001,1719PC0007,1719PC0009,1719PC0010,1719PC0015,1719PC0016		1478PC0006B,1478PC0008B,1478PC0013B,1478PC0014B,1478PC0015B,1478PC0016,1478PC0018,1478PC0020,1478PC0025,1719PC0001,1719PC0007,1719PC0009,1719PC0010,1719PC0015,1719PC0016" > exp
gemini query --header --show-samples -q "select chrom, start, end, ref, alt \
                                        from variants limit 5" test.query.db > obs
check obs exp
rm obs exp

##########################################################################
# 17. Test a query of the variants table with a
# request for all sample genotype types and a request for the sample names
##########################################################################
echo "    query.t17...\c"
echo "chrom	start	end	ref	alt	gt_types	variant_samples	HET_samples	HOM_ALT_samples
chr1	30547	30548	T	G	2,2,2,2,2,2,0,2,2,2,2,2,2,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,0,2,3,2,3,2,2,2,2,2,0,2,2,2,2,2,2,2	1478PC0016,1719PC0007,1719PC0009		1478PC0016,1719PC0007,1719PC0009
chr1	30859	30860	G	C	0,0,0,0,2,2,0,2,0,0,0,2,2,0,0,2,2,2,2,2,0,2,2,0,2,0,0,0,0,3,0,2,0,2,2,2,2,2,2,2,2,2,1,2,0,2,0,2,2,0,0,0,0,2,2,2,2,2,2,2	1478PC0017B,1719PC0005	1719PC0005	1478PC0017B
chr1	30866	30869	CCT	C	0,0,1,0,2,2,0,2,0,1,0,2,2,0,0,2,2,2,2,2,0,2,2,1,2,0,3,0,0,0,0,2,0,2,2,2,2,0,2,2,2,2,1,2,0,2,0,0,2,0,0,0,0,2,2,2,2,2,2,2	1094PC0012,1094PC0021,1478PC0011,1478PC0014B,1719PC0005	1094PC0012,1094PC0021,1478PC0011,1719PC0005	1478PC0014B
chr1	30894	30895	T	C	1,1,0,0,2,2,0,2,0,0,0,2,0,0,0,2,2,2,0,2,0,2,2,0,2,2,0,0,0,3,0,2,0,2,2,2,2,0,2,2,2,2,0,2,0,2,0,0,2,2,0,0,0,0,2,2,2,2,2,2	1094PC0005,1094PC0009,1478PC0017B	1094PC0005,1094PC0009	1478PC0017B
chr1	30922	30923	G	T	2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,3,2,2,2,2,3,3,3,3,2,3,2,3,2,2,2,2,3,3,2,2,2,2,2,3,2,3,3,2,2,2,2,3,3,2,2,2,2,2,2	1478PC0006B,1478PC0008B,1478PC0013B,1478PC0014B,1478PC0015B,1478PC0016,1478PC0018,1478PC0020,1478PC0025,1719PC0001,1719PC0007,1719PC0009,1719PC0010,1719PC0015,1719PC0016		1478PC0006B,1478PC0008B,1478PC0013B,1478PC0014B,1478PC0015B,1478PC0016,1478PC0018,1478PC0020,1478PC0025,1719PC0001,1719PC0007,1719PC0009,1719PC0010,1719PC0015,1719PC0016" > exp
gemini query --header --show-samples -q "select chrom, start, end, ref, alt, gt_types \
                                        from variants limit 5" test.query.db > obs
check obs exp
rm obs exp

####################################################################
# 18. Test tokenizing a SELECT list using spaces
####################################################################
echo "    query.t18...\c"
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
gemini query -q "select chrom, start, end, ref, alt from variants limit 10" test.query.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 19. Test tokenizing a SELECT list using spaces
####################################################################
echo "    query.t19...\c"
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
gemini query -q "select chrom,start,end,ref,alt from variants limit 10" test.query.db \
       > obs
check obs exp

####################################################################
# 20. Test tokenizing a SELECT list spaces and no spaces
####################################################################
echo "    query.t20...\c"
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
gemini query -q "select chrom, start,end, ref,alt from variants limit 10" test.query.db \
       > obs
check obs exp
rm obs exp

########################################################################
# 21. Test tokenizing a SELECT list spaces and no spaces and GT columns
########################################################################
echo "    query.t21...\c"
echo "chr1	30547	30548	T	G	T/T
chr1	30859	30860	G	C	G/G
chr1	30866	30869	CCT	C	CCT/CCT
chr1	30894	30895	T	C	T/T
chr1	30922	30923	G	T	./.
chr1	69269	69270	A	G	G/G
chr1	69427	69428	T	G	T/T
chr1	69510	69511	A	G	A/G
chr1	69760	69761	A	T	A/A
chr1	69870	69871	G	A	G/G" > exp
gemini query -q "select chrom, start,end, ref,alt,gts.1094PC0018 from variants limit 10" test.query.db \
       > obs
check obs exp
rm obs exp

########################################################################
# 22. Test transposed ped (TPED) query format
########################################################################
echo "    query.t22...\c"
echo "10 10:1142207-1142208:C 0 1142207 C C C C C C C C
10 10:48003991-48003992:C|T 0 48003991 T T C T C T C C
10 10:52004314-52004315:C 0 52004314 0 0 0 0 C C C C
10 10:52497528-52497529:C 0 52497528 0 0 C C C C 0 0
16 16:72057434-72057435:C|T 0 72057434 C T C C C C C C
10 10:126678091-126678092:A|G 0 126678091 G G G G G G G A
10 10:135210790-135210791:C|T 0 135210790 T T C C C C T T
10 10:135336655-135336656:A 0 135336655 0 0 A A 0 0 A A
10 10:135369531-135369532:C|T 0 135369531 T T T C T C T T" > exp
gemini query --format tped -q "select * from variants" test4.snpeff.ped.db > obs
check obs exp
rm obs exp

########################################################################
# 23. Test exclude-phenotype query
########################################################################
echo "    query.t23...\c"
echo "C/T,C/C,C/C,C/C	1,0,0,0
G/G,G/G,G/G,G/A	0,0,0,1" > exp
gemini query --sample-filter "phenotype!=2" --in only -q "select gts, gt_types from variants" test4.snpeff.ped.db > obs
#gemini query --exclude-phenotype affected -q "select gts, gt_types from variants" test4.snpeff.ped.db > obs
check obs exp
rm obs exp

########################################################################
# 24. Test phenotype query
########################################################################
echo "    query.t24...\c"
echo "./.,C/C,C/C,./.	2,3,3,2
T/T,C/C,C/C,T/T	0,3,3,0
T/T,T/C,T/C,T/T	0,1,1,0" > exp
gemini query --sample-filter "phenotype=2" --in only all -q "select gts, gt_types from variants" test4.snpeff.ped.db > obs
check obs exp
rm obs exp

########################################################################
# 25. Test region specific query
########################################################################
echo "    query.t25...\c"
echo "chr1	30859	G	30860	C" > exp
gemini query --region chr1:30859-30900 -q "select chrom, start, ref, end, alt from variants"  test1.snpeff.db > obs
check obs exp
rm obs exp

########################################################################
# 26. Test family-wise query
########################################################################
echo "    query.t26...\c"
echo "T/T,T/T,T/C,T/T,T/T,T/T,T/T,T/T,C/C	0,0,1,0,0,0,0,0,3
C/T,C/T,T/T,C/C,C/C,C/T,C/T,C/T,C/T	1,1,3,0,0,1,1,1,1
C/T,C/T,C/T,C/T,C/T,T/T,C/C,C/C,C/T	1,1,1,1,1,3,0,0,1
G/G,G/G,G/A,G/G,G/G,G/A,G/A,G/A,G/A	0,0,1,0,0,1,1,1,1
T/T,T/T,T/C,T/T,T/T,T/C,T/T,T/T,T/C	0,0,1,0,0,1,0,0,1" > exp
gemini query  --min-kindreds 2 --family-wise --sample-filter "phenotype=2" --in all -q "select gts, gt_types from variants" test.family.db > obs
check obs exp
rm obs exp

########################################################################
# 27. Test family-wise phenotype exclusion query
########################################################################
echo "    query.t27...\c"
echo "T/T,T/T,T/C,T/T,T/T,T/T,T/T,T/T,C/C	0,0,1,0,0,0,0,0,3
T/T,T/T,T/C,T/T,T/T,T/C,T/T,T/T,T/C	0,0,1,0,0,1,0,0,1
T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/C	0,0,0,0,0,0,0,0,1" > exp
gemini query  --in none --sample-filter "phenotype=1" -q "select gts, gt_types from variants" test.family.db > obs
check obs exp
rm obs exp

########################################################################
# 28. Test the extended ped sample-filter query
########################################################################
echo "    query.t28...\c"
echo "G/G,G/G,G/G,G/A	0,0,0,1" > exp
gemini query  --in only all --sample-filter "phenotype=1 and hair_color='blue'" -q "select gts, gt_types from variants" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 29. Test the carrier/noncarrier column summary
########################################################################
echo "    query.t29...\c"
echo "chrom	start	ref	alt	gt_types	variant_samples	HET_samples	HOM_ALT_samples	unaffected_carrier	affected_carrier	unaffected_noncarrier	affected_noncarrier	unknown
chr10	1142207	T	C	3,3,3,3	M10475,M10478,M10500,M128215		M10475,M10478,M10500,M128215	2	2	0	0	0
chr10	48003991	C	T	3,1,1,0	M10475,M10478,M10500	M10478,M10500	M10475	1	2	1	0	0
chr10	52004314	T	C	2,2,3,3	M10500,M128215		M10500,M128215	1	1	0	0	2
chr10	52497528	G	C	2,3,3,2	M10478,M10500		M10478,M10500	0	2	0	0	2
chr16	72057434	C	T	1,0,0,0	M10475	M10475		1	0	1	2	0
chr10	126678091	G	A	0,0,0,1	M128215	M128215		1	0	1	2	0
chr10	135210790	T	C	0,3,3,0	M10478,M10500		M10478,M10500	0	2	2	0	0
chr10	135336655	G	A	2,3,2,3	M10478,M128215		M10478,M128215	1	1	0	0	2
chr10	135369531	T	C	0,1,1,0	M10478,M10500	M10478,M10500		0	2	2	0	0" > exp
gemini query --show-samples --carrier-summary-by-phenotype affected --header -q "select chrom, start, ref, alt, gt_types from variants" extended_ped.db > obs
check obs exp
rm obs exp

#########################################################################
# 30. Test a basic query of the variants table with show-families
#########################################################################
echo "    query.t30...\c"
echo "chrom	start	end	ref	alt	families
chr10	1142207	1142208	T	C	1
chr10	48003991	48003992	C	T	1
chr10	52004314	52004315	T	C	1
chr10	52497528	52497529	G	C	1
chr16	72057434	72057435	C	T	1" > exp
gemini query --header --show-families -q "select chrom, start, end, ref, alt from variants limit 5" extended_ped.db > obs
check obs exp
rm obs exp

####################################################################
# 31. Test that rows are filtered based on a --gt-filter if
#     a GT* column is SELECTed.
####################################################################
echo "    query.t31...\c"
echo "chr1	69269	69270	A	G	OR4F5	3
chr1	69510	69511	A	G	OR4F5	3
chr1	69760	69761	A	T	OR4F5	3
chr1	861629	861630	G	A	SAMD11	3
chr1	861807	861808	A	G	SAMD11	3
chr1	866318	866319	G	A	SAMD11	3
chr1	866510	866511	C	CCCCT	SAMD11	3
chr1	866892	866893	T	C	SAMD11	3
chr1	866919	866920	A	G	SAMD11	3
chr1	870902	870903	T	C	SAMD11	3" > exp
gemini query -q "select chrom, start, end, ref, alt, gene, gt_types.1094PC0019 \
                 from variants" \
             --gt-filter "gt_types.1094PC0019 == HOM_ALT" test.query.db | head \
       > obs
check obs exp
rm obs exp

####################################################################
# 32. Test that rows are filtered based on a --gt-filter if
#     a GT* column is NOT SELECTed.
####################################################################
echo "    query.t32...\c"
echo "chr1	69269	69270	A	G	OR4F5
chr1	69510	69511	A	G	OR4F5
chr1	69760	69761	A	T	OR4F5
chr1	861629	861630	G	A	SAMD11
chr1	861807	861808	A	G	SAMD11
chr1	866318	866319	G	A	SAMD11
chr1	866510	866511	C	CCCCT	SAMD11
chr1	866892	866893	T	C	SAMD11
chr1	866919	866920	A	G	SAMD11
chr1	870902	870903	T	C	SAMD11" > exp
gemini query -q "select chrom, start, end, ref, alt, gene \
                 from variants" \
             --gt-filter "gt_types.1094PC0019 == HOM_ALT" test.query.db | head \
       > obs
check obs exp
rm obs exp

####################################################################
# 33. Test that non-genotype columns that contain the substring "gt"
# execute properly
####################################################################
echo "    query.t33...\c"
echo "chrom	start	end	ref	alt	aa_length	gene
chr1	30547	30548	T	G	85	FAM138A
chr1	30859	30860	G	C	85	FAM138A
chr1	30866	30869	CCT	C	85	FAM138A
chr1	30894	30895	T	C	85	FAM138A
chr1	30922	30923	G	T	85	FAM138A
chr1	69269	69270	A	G	305	OR4F5
chr1	69427	69428	T	G	305	OR4F5
chr1	69510	69511	A	G	305	OR4F5
chr1	69760	69761	A	T	305	OR4F5" > exp
gemini query --header -q "select chrom, start, end, ref, alt, aa_length, gene \
                 from variants" test.query.db | head \
       > obs
check obs exp
rm obs exp