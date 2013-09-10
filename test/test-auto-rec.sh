###################################################################
# 1. Test basic auto_recessive functionality
###################################################################
echo "    auto_rec.t1...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/T,C/T,T/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	C/T,C/T,T/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/C,T/C,C/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/C,T/C,C/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/C,T/C,C/C	59,49,64	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    test.auto_rec.db > obs
check obs exp
rm obs exp


###################################################################
# 2. Test with a minimum of 2 kindreds
###################################################################
echo "    auto_rec.t2...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/T,C/T,T/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	C/T,C/T,T/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/C,T/C,C/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/C,T/C,C/C	59,49,64	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_rec.db > obs
check obs exp
rm obs exp


###################################################################
# 3. Test with a minimum of 3 kindreds
###################################################################
echo "    auto_rec.t3...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    test.auto_rec.db > obs
check obs exp
rm obs exp


###################################################################
# 4. Test with a filter and minimum of 2 kindreds, HIGH severity
###################################################################
echo "    auto_rec.t4...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/C,T/C,C/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/C,T/C,C/C	59,49,64	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    --filter "impact_severity = 'HIGH'" \
    test.auto_rec.db > obs
check obs exp
rm obs exp


###################################################################
# 5. Test with a filter and minimum of 1 kindreds, HIGH severity and min depth of 40
###################################################################
echo "    auto_rec.t5...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/C,T/C,C/C	59,49,64	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    --filter "impact_severity = 'HIGH'" \
    -d 40 \
    test.auto_rec.db > obs
check obs exp
rm obs exp
