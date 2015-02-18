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
# 1. Test basic auto_dominant functionality
###################################################################
echo "    auto_dom.t1...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	G/A,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    test.auto_dom.db > obs
check obs exp
rm obs exp


###################################################################
# 2. Test with a minimum of 2 kindreds
###################################################################
echo "    auto_dom.t2...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
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
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
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
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    --filter "impact_severity = 'HIGH'" \
    test.auto_dom.db > obs
check obs exp
rm obs exp

###################################################################
# 5. Test with a filter and minimum of 1 kindreds, min depth 40
###################################################################
echo "    auto_dom.t5...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    --filter "impact_severity = 'HIGH'" \
    -d 40 \
    test.auto_dom.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test with one family lacking parents, requiring one kindred.
###################################################################
echo "    auto_dom.t6...\c"
echo "WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	G/A,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_dom.no_parents.db &> obs
check obs exp
#rm obs exp

####################################################################
## 7. Test with one family lacking parents, requiring two kindreds.
####################################################################
echo "    auto_dom.t7...\c"
echo "WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_dom.no_parents.db &> obs
check obs exp
rm obs exp

###################################################################
# 8. Test with two families lacking parents, requiring one kindred.
###################################################################
echo "    auto_dom.t8...\c"
echo "WARNING: Unable to identify parents for family (3). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	G/A,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(affected),3_mom(unknown),3_kid(affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_dom.no_parents.2.db &> obs
check obs exp
rm obs exp

###################################################################
# 9. Test with two families lacking parents, requiring two kindreds.
###################################################################
echo "    auto_dom.t9...\c"
echo "WARNING: Unable to identify parents for family (3). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(affected),3_mom(unknown),3_kid(affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_dom.no_parents.2.db &> obs
check obs exp
rm obs exp

###################################################################
# 10. Test with three family lacking parents, requiring one kindred.
###################################################################
echo "    auto_dom.t10...\c"
echo "WARNING: Unable to identify parents for family (1). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (3). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	G/A,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	T/T,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(affected),3_mom(unknown),3_kid(affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_dom.no_parents.3.db &> obs
check obs exp
rm obs exp

###################################################################
# 11. Test with three family lacking parents, requiring two kindreds.
###################################################################
echo "    auto_dom.t11...\c"
echo "WARNING: Unable to identify parents for family (1). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (3). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	T/T,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(affected),3_mom(unknown),3_kid(affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_dom.no_parents.3.db &> obs
check obs exp
rm obs exp

###################################################################
# 12. Test with three family lacking parents, requiring two kindreds.
###################################################################
echo "    auto_dom.t12...\c"
echo "WARNING: Unable to identify parents for family (1). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (3). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
WARNING: Unable to identify parents for family (2). Consequently, GEMINI will solely place genotype requirements on subjects based on their phenotype.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(affected),3_mom(unknown),3_kid(affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	C/C,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
1	1_dad(unaffected),1_mom(unaffected),1_kid(affected)	T/T,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
2	2_dad(unaffected),2_mom(affected),2_kid(affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(affected),3_mom(unknown),3_kid(affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    test.auto_dom.no_parents.3.db &> obs
check obs exp
rm obs exp


###################################################################
# 13. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_dom.t13...\c"
echo "WARNING: Unable to identify at least one affected individual for family (1). Consequently, GEMINI will not screen for variants in this family.
WARNING: Unable to identify at least one affected individual for family (3). Consequently, GEMINI will not screen for variants in this family.
WARNING: Unable to identify at least one affected individual for family (2). Consequently, GEMINI will not screen for variants in this family.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_dom.no_parents.4.db &> obs
check obs exp
rm obs exp

###################################################################
# 14. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_dom.t14...\c"
echo "WARNING: Unable to identify at least one affected individual for family (1). Consequently, GEMINI will not screen for variants in this family.
WARNING: Unable to identify at least one affected individual for family (3). Consequently, GEMINI will not screen for variants in this family.
WARNING: Unable to identify at least one affected individual for family (2). Consequently, GEMINI will not screen for variants in this family.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_dom.no_parents.4.db &> obs
check obs exp
rm obs exp

###################################################################
# 15. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_dom.t15...\c"
echo "WARNING: Unable to identify at least one affected individual for family (1). Consequently, GEMINI will not screen for variants in this family.
WARNING: Unable to identify at least one affected individual for family (3). Consequently, GEMINI will not screen for variants in this family.
WARNING: Unable to identify at least one affected individual for family (2). Consequently, GEMINI will not screen for variants in this family.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    test.auto_dom.no_parents.4.db &> obs
check obs exp
rm obs exp

###################################################################
# 16. Test with unknown parents, yet affected kid.
###################################################################
echo "    auto_dom.t16...\c"
echo "WARNING: Unable to identify at least one affected individual for family (3). Consequently, GEMINI will not screen for variants in this family.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
2	2_dad(father; unknown),2_mom(mother; unknown),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unknown),2_mom(mother; unknown),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unknown),2_mom(mother; unknown),2_kid(child; affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_dom.no_parents.5.db &> obs
check obs exp
rm obs exp

###################################################################
# 17. Test with --families
###################################################################
echo "    auto_dom.t17...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	C/C,C/T,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	G/A,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
2	2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)	T/T,T/C,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --families 2,3 \
    test.auto_dom.db &> obs
check obs exp
rm obs exp

###################################################################
# 18. Test with --families
###################################################################
echo "    auto_dom.t18...\c"
echo "family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	C/T,C/C,C/T	39,29,24	ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	G/A,G/G,G/A	39,29,24	SPRN	chr10	135336655	135336656	G	A	intron	LOW
3	3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)	T/C,T/T,T/C	39,29,24	WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --families 3 \
    test.auto_dom.db &> obs
check obs exp
rm obs exp

###################################################################
# 19. Test with --families but not valid
###################################################################
echo "    auto_dom.t19...\c"
echo "WARNING: Unable to identify at least one affected individual for family (3). Consequently, GEMINI will not screen for variants in this family.
family_id	family_members	family_genotypes	family_genotype_depths	gene	chrom	start	end	ref	alt	impact	impact_severity" > exp
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --families 3 \
    test.auto_dom.no_parents.5.db &> obs
check obs exp
rm obs exp