
echo "test-dashes.t01"
gemini query -q "select chrom,start,end,gt_types.Na-12877 from variants" test.dashes.db > obs

echo "chr1	10670	10671	1
chr1	28493	28494	3
chr1	28627	28628	1
chr1	137122	137123	0
chr1	267558	267560	1
chr1	537540	537541	1
chr1	537969	537970	1
chr1	540709	540712	1
chr1	541051	541052	1
chr1	547422	547423	1
chr1	547518	547519	1
chr1	547752	547753	1
chr1	548245	548246	1
chr1	589081	589086	3
chr1	589081	589136	2" > exp
check obs exp



echo "test-dashes.t02"
gemini mendel_errors --columns "end,gt_types.Na-12877" test.dashes.db > obs
echo "end	variant_id	gt_types.Na-12877	family_id	family_members	family_genotypes	samples	family_count	violation	violation_prob
10671	1	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	G/G,G/G,G/C	Na-12877	1	plausible de novo	0.962
28494	2	3	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	T/C,T/T,C/C	Na-12877	1	loss of heterozygosity	0.660
28628	3	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	C/C,C/C,C/T	Na-12877	1	plausible de novo	0.989
137123	4	0	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	T/T,G/G,T/T	Na-12877	1	uniparental disomy	0.881
267560	5	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	C/C,C/C,CT/C	Na-12877	1	plausible de novo	0.896
537541	6	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	T/T,T/T,T/TC	Na-12877	1	plausible de novo	0.729
537970	7	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	C/C,C/C,C/T	Na-12877	1	plausible de novo	0.928
540712	8	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	CGA/CGA,CGA/CGA,CGA/C	Na-12877	1	plausible de novo	0.814
541052	9	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	T/T,T/T,T/C	Na-12877	1	plausible de novo	0.750
547423	10	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	G/G,G/G,G/A	Na-12877	1	plausible de novo	0.746
547519	11	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	G/G,G/G,G/T	Na-12877	1	plausible de novo	1.000
547753	12	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	G/G,G/G,G/A	Na-12877	1	plausible de novo	0.754
548246	13	1	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	T/T,T/T,T/TA	Na-12877	1	plausible de novo	0.925
589086	14	3	CEPH1463	Na-12889(dad;unknown),Na-12890(mom;unknown),Na-12877(child;unknown)	G/G,GAGAA/GAGAA,G/G	Na-12877	1	uniparental disomy	0.940" > exp
check obs exp


echo "test-dashes.t03"
gemini query -q "select chrom,start,end,gt_types.Na-12877 from variants" test.dashes.db --gt-filter "gt_types.Na-12877 == HOM_REF" > obs
echo "chr1	137122	137123	0" > exp
check obs exp

gemini bcolz_index test.dashes.db

echo "test-dashes.t04"
gemini query -q "select chrom,start,end,gt_types.Na-12877 from variants" test.dashes.db --gt-filter "gt_types.Na-12877 == HOM_REF" --use-bcolz > obs
echo "chr1	137122	137123	0" > exp
check obs exp

echo "test-dashes.t05"
gemini comp_hets --column "chrom,start,gt_types.Na-12877" test.dashes.db > obs
echo "chrom	start	variant_id	gt_types.Na-12877	family_id	family_members	family_genotypes	samples	family_count	comp_het_id" > exp
check obs exp


echo "test-dashes.t06"
gemini autosomal_recessive --column "gt_types.Na-12877" test.dashes.db > obs
echo "gene	variant_id	gt_types.Na-12877	family_id	family_members	family_genotypes	samples	family_count" > exp
