###################################################################
# 1. Test basic de_novo functionality
###################################################################
echo "    de_novo.t1...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    test.de_novo.db > obs
check obs exp
rm obs exp


###################################################################
# 2. Test de_novo with filter
###################################################################
echo "    de_novo.t2...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --filter "impact_severity = 'HIGH'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test de_novo with filter and minimum depth requirement
###################################################################
echo "    de_novo.t3...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --filter "impact_severity = 'MED'" \
    -d 40 \
    test.de_novo.db > obs
check obs exp
rm obs exp

