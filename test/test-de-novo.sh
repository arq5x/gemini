check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}
export -f check

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

###################################################################
# 4. Test de_novo with filter with min-kindreds
###################################################################
echo "    de_novo.t4...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 5. Test de_novo with filter with min-kindreds 2 and filter
###################################################################
echo "    de_novo.t5...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    --filter "impact_severity != 'LOW'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test de_novo with filter with min-kindreds 3 and filter
###################################################################
echo "    de_novo.t6...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    --filter "impact_severity != 'LOW'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 7. Test de_novo with filter with min-kindreds 1 and filter
###################################################################
echo "    de_novo.t7...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)	T/T,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    --filter "impact_severity != 'LOW'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 8. Test de_novo without --only-affected
###################################################################
echo "    de_novo.t8...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	T/T,T/T,T/C,T/C	39,29,24,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	G/G,G/G,G/A,G/G	39,29,24,39	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	T/T,T/T,T/T,T/C	39,29,24,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp

###################################################################
# 9. Test de_novo with --only-affected
###################################################################
echo "    de_novo.t9...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	T/T,T/T,T/C,T/C	39,29,24,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	T/T,T/T,T/T,T/C	39,29,24,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
3	3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)	T/T,T/T,T/C	50,50,50	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --only-affected \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp

###################################################################
# 10. Test de_novo with --only-affected and --families
###################################################################
echo "    de_novo.t10...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	T/T,T/T,T/C,T/C	39,29,24,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
1	1_dad(father; unaffected),1_mom(mother; unaffected),1_kid1(child; unaffected),1_kid2(child; affected)	T/T,T/T,T/T,T/C	39,29,24,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --only-affected \
    --families 1,2 \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp

###################################################################
# 11. Test de_novo with --only-affected and --families
###################################################################
echo "    de_novo.t11...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	G/G,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)	T/T,T/T,T/C	39,29,24	SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED" > exp
gemini de_novo  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --only-affected \
    --families 2 \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp
