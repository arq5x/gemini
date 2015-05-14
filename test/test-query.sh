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
echo "1	0	1094PC0005	0	0	-9	-9
2	0	1094PC0009	0	0	-9	-9
3	0	1094PC0012	0	0	-9	-9
4	0	1094PC0013	0	0	-9	-9
5	0	1094PC0016	0	0	-9	-9
6	0	1094PC0017	0	0	-9	-9
7	0	1094PC0018	0	0	-9	-9
8	0	1094PC0019	0	0	-9	-9
9	0	1094PC0020	0	0	-9	-9
10	0	1094PC0021	0	0	-9	-9" > exp
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
echo "chr1	30547	30548	None	1	1	T	G	50.09	None	snp	tv	0.116666666667
chr1	30859	30860	None	2	1	G	C	54.3	None	snp	tv	0.433333333333
chr1	30866	30869	None	3	1	CCT	C	49.48	None	indel	del	0.466666666667
chr1	30894	30895	None	4	1	T	C	51.79	None	snp	ts	0.483333333333
chr1	30922	30923	None	5	1	G	T	601.49	None	snp	tv	0.25" > exp
gemini query -q "select * \
                 from variants \
                 limit 5" test.query.db | cut -f1-13 \
       > obs
check obs exp
rm obs exp


##############################################################################
# 11. Test a query of the variants table with a select * and a genotype column
###############################################################################
echo "    query.t11...\c"
echo "chr1	30547	30548	T	G	T/T
chr1	30859	30860	G	C	G/G
chr1	30866	30869	CCT	C	CCT/CCT
chr1	30894	30895	T	C	T/T
chr1	30922	30923	G	T	./." > exp
gemini query -q "select *, gts.1094PC0018 \
                  from variants \
                  limit 5" test.query.db | awk '{OFS="\t"}{print $1,$2,$3,$7,$8,$NF}' > obs
check obs exp
rm obs exp


##############################################################################
# 12. Test a query of the variants table with a select * and the full genotype column
###############################################################################
echo "    query.t12...\c"
echo "chr1	30547	30548	T	G	./.,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,G/G,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,T/T,./.,G/G,./.,G/G,./.,./.,./.,./.,./.,T/T,./.,./.,./.,./.,./.,./.,./.
chr1	30859	30860	G	C	G/G,G/G,G/G,G/G,./.,./.,G/G,./.,G/G,G/G,G/G,./.,./.,G/G,G/G,./.,./.,./.,./.,./.,G/G,./.,./.,G/G,./.,G/G,G/G,G/G,G/G,C/C,G/G,./.,G/G,./.,./.,./.,./.,./.,./.,./.,./.,./.,G/C,./.,G/G,./.,G/G,./.,./.,G/G,G/G,G/G,G/G,./.,./.,./.,./.,./.,./.,./.
chr1	30866	30869	CCT	C	CCT/CCT,CCT/CCT,CCT/C,CCT/CCT,./.,./.,CCT/CCT,./.,CCT/CCT,CCT/C,CCT/CCT,./.,./.,CCT/CCT,CCT/CCT,./.,./.,./.,./.,./.,CCT/CCT,./.,./.,CCT/C,./.,CCT/CCT,C/C,CCT/CCT,CCT/CCT,CCT/CCT,CCT/CCT,./.,CCT/CCT,./.,./.,./.,./.,CCT/CCT,./.,./.,./.,./.,CCT/C,./.,CCT/CCT,./.,CCT/CCT,CCT/CCT,./.,CCT/CCT,CCT/CCT,CCT/CCT,CCT/CCT,./.,./.,./.,./.,./.,./.,./.
chr1	30894	30895	T	C	T/C,T/C,T/T,T/T,./.,./.,T/T,./.,T/T,T/T,T/T,./.,T/T,T/T,T/T,./.,./.,./.,T/T,./.,T/T,./.,./.,T/T,./.,./.,T/T,T/T,T/T,C/C,T/T,./.,T/T,./.,./.,./.,./.,T/T,./.,./.,./.,./.,T/T,./.,T/T,./.,T/T,T/T,./.,./.,T/T,T/T,T/T,T/T,./.,./.,./.,./.,./.,./.
chr1	30922	30923	G	T	./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,./.,T/T,./.,T/T,./.,./.,./.,./.,T/T,T/T,T/T,T/T,./.,T/T,./.,T/T,./.,./.,./.,./.,T/T,T/T,./.,./.,./.,./.,./.,T/T,./.,T/T,T/T,./.,./.,./.,./.,T/T,T/T,./.,./.,./.,./.,./.,./." > exp
gemini query -q "select *, gts \
                  from variants \
                  limit 5" test.query.db | awk '{OFS="\t"}{print $1,$2,$3,$7,$8,$NF}' > obs
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
# NOTE: this test fails with --use-bcolz because of the the limit 5
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
echo "chrom	start	end	ref	alt	variant_samples	het_samples	hom_alt_samples
chr1	30547	30548	T	G	1478PC0016,1719PC0007,1719PC0009		1478PC0016,1719PC0007,1719PC0009
chr1	30859	30860	G	C	1719PC0005,1478PC0017B	1719PC0005	1478PC0017B
chr1	30866	30869	CCT	C	1094PC0012,1094PC0021,1478PC0011,1719PC0005,1478PC0014B	1094PC0012,1094PC0021,1478PC0011,1719PC0005	1478PC0014B
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
echo "chrom	start	end	ref	alt	gt_types	variant_samples	het_samples	hom_alt_samples
chr1	30547	30548	T	G	2,2,2,2,2,2,0,2,2,2,2,2,2,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,0,2,3,2,3,2,2,2,2,2,0,2,2,2,2,2,2,2	1478PC0016,1719PC0007,1719PC0009		1478PC0016,1719PC0007,1719PC0009
chr1	30859	30860	G	C	0,0,0,0,2,2,0,2,0,0,0,2,2,0,0,2,2,2,2,2,0,2,2,0,2,0,0,0,0,3,0,2,0,2,2,2,2,2,2,2,2,2,1,2,0,2,0,2,2,0,0,0,0,2,2,2,2,2,2,2	1719PC0005,1478PC0017B	1719PC0005	1478PC0017B
chr1	30866	30869	CCT	C	0,0,1,0,2,2,0,2,0,1,0,2,2,0,0,2,2,2,2,2,0,2,2,1,2,0,3,0,0,0,0,2,0,2,2,2,2,0,2,2,2,2,1,2,0,2,0,0,2,0,0,0,0,2,2,2,2,2,2,2	1094PC0012,1094PC0021,1478PC0011,1719PC0005,1478PC0014B	1094PC0012,1094PC0021,1478PC0011,1719PC0005	1478PC0014B
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
echo "10 1 0 1142207 C C C C C C C C
10 2 0 48003991 T T C T C T C C
10 3 0 52004314 0 0 0 0 C C C C
10 4 0 52497528 0 0 C C C C 0 0
10 5 0 126678091 G G G G G G G A
10 6 0 135210790 T T C C C C T T
10 7 0 135336655 0 0 A A 0 0 A A
10 8 0 135369531 T T T C T C T T
16 9 0 72057434 C T C C C C C C" > exp
gemini query --format tped -q "select * from variants" test4.snpeff.ped.db > obs
check obs exp
rm obs exp

########################################################################
# 23. Test exclude-phenotype query
########################################################################
echo "    query.t23...\c"
echo "G/G,G/G,G/G,G/A	0,0,0,1
C/T,C/C,C/C,C/C	1,0,0,0" > exp
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
echo "chrom	start	ref	alt	gt_types	variant_samples	het_samples	hom_alt_samples	unaffected_carrier	affected_carrier	unaffected_noncarrier	affected_noncarrier	unknown
chr10	1142207	T	C	3,3,3,3	M10475,M10478,M10500,M128215		M10475,M10478,M10500,M128215	2	2	0	0	0
chr10	48003991	C	T	3,1,1,0	M10478,M10500,M10475	M10478,M10500	M10475	1	2	1	0	0
chr10	52004314	T	C	2,2,3,3	M10500,M128215		M10500,M128215	1	1	0	0	2
chr10	52497528	G	C	2,3,3,2	M10478,M10500		M10478,M10500	0	2	0	0	2
chr10	126678091	G	A	0,0,0,1	M128215	M128215		1	0	1	2	0
chr10	135210790	T	C	0,3,3,0	M10478,M10500		M10478,M10500	0	2	2	0	0
chr10	135336655	G	A	2,3,2,3	M10478,M128215		M10478,M128215	1	1	0	0	2
chr10	135369531	T	C	0,1,1,0	M10478,M10500	M10478,M10500		0	2	2	0	0
chr16	72057434	C	T	1,0,0,0	M10475	M10475		1	0	1	2	0" > exp
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
chr10	126678091	126678092	G	A	1" > exp
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
                 from variants limit 9" test.query.db \
       > obs
check obs exp
rm obs exp

#########################################################################
# 34. Test the gene_detailed table and the join on variants table
#########################################################################
echo "    query.t34...\c"
echo "variant_id	chrom	gene	transcript_status	transcript	transcript_start	transcript_end	synonym	rvis_pct	protein_length	impact
46	chr1	SAMD11	KNOWN	ENST00000342066	861118	879955	MGC45873	None	681	frame_shift
578	chr1	TNFRSF18	PUTATIVE	ENST00000486728	1139224	1141060	AITR,CD357,GITR	None	169	frame_shift
733	chr1	SCNN1D	NOVEL	ENST00000470022	1217305	1221548	ENaCdelta,dNaCh	96.77990092	138	stop_gain" > exp

gemini query --header -q "select v.variant_id, v.chrom, v.gene, \
	           g.transcript_status, g.transcript, g.transcript_start, \
			    g.transcript_end, g.synonym, g.rvis_pct, g.protein_length, \
				v.impact from variants v, gene_detailed g \
					
					WHERE v.chrom = g.chrom AND \
						  v.gene = g.gene AND v.impact_severity='HIGH' AND \
						  v.biotype='protein_coding' AND \
						  v.transcript = g.transcript" test.query.db > obs
check obs exp
rm obs exp

#########################################################################
# 35. Test the gene_detailed table and the join on variant_impacts table
#########################################################################
echo "    query.t35...\c"
echo "gene	transcript_status	transcript	transcript_start	transcript_end	synonym	rvis_pct	protein_length	impact
SAMD11	KNOWN	ENST00000342066	861118	879955	MGC45873	None	681	frame_shift
TNFRSF18	PUTATIVE	ENST00000486728	1139224	1141060	AITR,CD357,GITR	None	169	frame_shift
TNFRSF18	KNOWN	ENST00000379265	1139224	1141951	AITR,CD357,GITR	None	234	frame_shift
TNFRSF18	KNOWN	ENST00000379268	1138891	1142071	AITR,CD357,GITR	None	241	frame_shift
TNFRSF18	KNOWN	ENST00000328596	1138888	1141951	AITR,CD357,GITR	None	255	frame_shift
SCNN1D	NOVEL	ENST00000470022	1217305	1221548	ENaCdelta,dNaCh	96.77990092	138	stop_gain
SCNN1D	NOVEL	ENST00000470022	1217305	1221548	ENaCdelta,dNaCh	96.77990092	138	frame_shift
SCNN1D	KNOWN	ENST00000325425	1217489	1227404	ENaCdelta,dNaCh	96.77990092	704	frame_shift
SCNN1D	KNOWN	ENST00000379116	1215816	1227399	ENaCdelta,dNaCh	96.77990092	802	frame_shift
SCNN1D	KNOWN	ENST00000338555	1215968	1227404	ENaCdelta,dNaCh	96.77990092	638	frame_shift
SCNN1D	KNOWN	ENST00000400928	1217576	1227409	ENaCdelta,dNaCh	96.77990092	638	frame_shift" > exp

gemini query --header -q "select v.gene, g.transcript_status,g.transcript, g.transcript_start, \
	g.transcript_end, g.synonym, g.rvis_pct, g.protein_length, \
    v.impact from variant_impacts v, gene_detailed g \
		WHERE v.transcript = g.transcript AND \
              v.gene = g.gene AND \
	          v.impact_severity='HIGH' AND \
              v.biotype='protein_coding'" test.query.db > obs

check obs exp
rm obs exp

###########################################################################
# 36. Test the gene_summary table and the join on variants table
###########################################################################
echo "    query.t36...\c"
echo "chrom	gene	strand	transcript_min_start	transcript_max_end	synonym	rvis_pct	impact
chr1	SAMD11	1	860260	879955	MGC45873	None	frame_shift
chr1	TNFRSF18	-1	1138888	1142071	AITR,CD357,GITR	None	frame_shift
chr1	SCNN1D	1	1215816	1227409	ENaCdelta,dNaCh	96.77990092	stop_gain
chr1	SCNN1D	1	1215816	1227409	ENaCdelta,dNaCh	96.77990092	frame_shift" > exp

gemini query --header -q "select v.chrom, v.gene, g.strand, g.transcript_min_start, g.transcript_max_end, \
g.synonym, g.rvis_pct, v.impact from variants v, gene_summary g \
WHERE v.chrom = g.chrom AND \
v.gene = g.gene AND \
v.impact_severity='HIGH'" test.query.db > obs

check obs exp
rm obs exp 

############################################################################
# 37. Test the gene_summary table and the join on variant_impacts table
############################################################################
echo "    query.t37...\c"
echo "gene	impact	transcript	transcript_min_start	transcript_max_end	rvis_pct	synonym
SCNN1D	stop_gain	ENST00000470022	1215816	1227409	96.77990092	ENaCdelta,dNaCh" > exp

gemini query --header -q "select g.gene, v.impact, v.transcript, \
	   g.transcript_min_start, g.transcript_max_end, g.rvis_pct, g.synonym \
		  from gene_summary g, variant_impacts v \
			  where g.gene=v.gene AND \
				g.gene ='SCNN1D' AND \
				v.impact ='stop_gain'" test.query.db > obs
check obs exp
rm obs exp


############################################################################
# 38. Test the loading of gene_summary table with multiple cores
############################################################################
echo "    query.t38...\c"
echo "chrom	gene	strand	transcript_min_start	transcript_max_end	synonym	rvis_pct	impact
chr1	SAMD11	1	860260	879955	MGC45873	None	frame_shift
chr1	TNFRSF18	-1	1138888	1142071	AITR,CD357,GITR	None	frame_shift
chr1	SCNN1D	1	1215816	1227409	ENaCdelta,dNaCh	96.77990092	stop_gain
chr1	SCNN1D	1	1215816	1227409	ENaCdelta,dNaCh	96.77990092	frame_shift" > exp

gemini query --header -q "select v.chrom, v.gene, g.strand, g.transcript_min_start, g.transcript_max_end, \
g.synonym, g.rvis_pct, v.impact from variants v, gene_summary g \
WHERE v.chrom = g.chrom AND \
v.gene = g.gene AND \
v.impact_severity='HIGH'" test.query.core.db > obs

check obs exp
rm obs exp 

###########################################################################
# 39. Test the loading of gene_detailed table with multiple cores
###########################################################################
echo "    query.t39...\c"
echo "variant_id	chrom	gene	transcript_status	transcript	transcript_start	transcript_end	synonym	rvis_pct	protein_length	impact
46	chr1	SAMD11	KNOWN	ENST00000342066	861118	879955	MGC45873	None	681	frame_shift
578	chr1	TNFRSF18	PUTATIVE	ENST00000486728	1139224	1141060	AITR,CD357,GITR	None	169	frame_shift
733	chr1	SCNN1D	NOVEL	ENST00000470022	1217305	1221548	ENaCdelta,dNaCh	96.77990092	138	stop_gain" > exp

gemini query --header -q "select v.variant_id, v.chrom, v.gene, \
	           g.transcript_status, g.transcript, g.transcript_start, \
			    g.transcript_end, g.synonym, g.rvis_pct, g.protein_length, \
				v.impact from variants v, gene_detailed g \
					
					WHERE v.chrom = g.chrom AND \
						  v.gene = g.gene AND v.impact_severity='HIGH' AND \
						  v.biotype='protein_coding' AND \
						  v.transcript = g.transcript" test.query.core.db > obs
check obs exp
rm obs exp

############################################################################
# 40. Test phenotype column of the gene_table
############################################################################
echo "    query.t40...\c"
echo "chrom	end	gene	mam_phenotype_id
chr1	949422	ISG15	MP:0005390
chr1	949608	ISG15	MP:0005390
chr1	949655	ISG15	MP:0005390
chr1	949832	ISG15	MP:0005390
chr1	977356	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	978762	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	979594	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	979748	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	981131	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	982205	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	983005	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	984711	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	984957	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	985378	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	985826	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	986885	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	990213	AGRN	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
chr1	1139216	TNFRSF18	MP:0005397,MP:0005384,MP:0005387
chr1	1139498	TNFRSF18	MP:0005397,MP:0005384,MP:0005387
chr1	1140813	TNFRSF18	MP:0005397,MP:0005384,MP:0005387
chr1	1147337	TNFRSF4	MP:0005384,MP:0005397,MP:0005378,MP:0002873,MP:0005388,MP:0005370,MP:0005387,MP:0005381
chr1	1149480	TNFRSF4	MP:0005384,MP:0005397,MP:0005378,MP:0002873,MP:0005388,MP:0005370,MP:0005387,MP:0005381" > exp

gemini query --header -q "select v.chrom, v.end, v.gene, g.mam_phenotype_id from variants v, \
	                      gene_summary g where v.chrom=g.chrom and v.gene=g.gene and \
						  v.impact_severity !='LOW' and mam_phenotype_id !='None'" test.query.db > obs

check obs exp
rm obs exp

#########################################################################
# 41. Show an expanded version of sample information with --format sampledetail
#########################################################################
echo "    query.t41...\c"
echo "chrom	start	ref	family_id	name	paternal_id	maternal_id	sex	phenotype
chr1	30547	T	0	1478PC0016	0	0	-9	-9
chr1	30547	T	0	1719PC0007	0	0	-9	-9
chr1	30547	T	0	1719PC0009	0	0	-9	-9" > exp
gemini query --header --format sampledetail --show-samples -q "select chrom, start, ref \
                                                                from variants limit 1" test.query.db > obs
check obs exp
rm obs exp










 
