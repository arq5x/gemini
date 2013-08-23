###################################################################
# 1. Test basic auto_dominant functionality
###################################################################
echo "    auto_dom.t1...\c"
echo "family_id	family_members	family_genotypes	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	G/A,G/G,G/A	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    test.auto_dom.db > obs
check obs exp
rm obs exp


###################################################################
# 2. Test with a minimum of 2 kindreds
###################################################################
echo "    auto_dom.t2...\c"
echo "family_id	family_members	family_genotypes	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_dom.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test with a minimum of 3 kindreds
###################################################################
echo "    auto_dom.t3...\c"
echo "family_id	family_members	family_genotypes	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    test.auto_dom.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test with a filter and minimum of 2 kindreds
###################################################################
echo "    auto_dom.t4...\c"
echo "family_id	family_members	family_genotypes	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    --filter "impact_severity = 'HIGH'" \
    test.auto_dom.db > obs
check obs exp
rm obs exp
