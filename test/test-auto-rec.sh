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
# 1. Test basic auto_recessive functionality
###################################################################
echo "    auto_rec.t1...\c"
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/C,T/C,C/C	3_kid	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp

gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    test.auto_rec.db > obs
check obs exp

rm obs exp

###################################################################
# 2. Test with a minimum of 2 kindreds
###################################################################
echo "    auto_rec.t2...\c"
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
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
touch exp
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
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
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
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	1" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    --filter "impact_severity = 'HIGH'" \
    -d 40 \
    test.auto_rec.db > obs
check obs exp
rm obs exp


###################################################################
# 6. Test with one family lacking parents, requiring one kindred.
###################################################################
echo "    auto_rec.t6...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/C,T/C,C/C	3_kid	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_rec.no_parents.db &> obs
check obs exp
rm obs exp

###################################################################
# 7. Test with one family lacking parents, requiring two kindreds.
###################################################################
echo "    auto_rec.t7...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_rec.no_parents.db &> obs
check obs exp
rm obs exp

###################################################################
# 8. Test with two families lacking parents, requiring one kindred.
###################################################################
echo "    auto_rec.t8...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
WARNING: auto-recessive called on family 2 where no affected has parents
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/C,T/C,C/C	3_kid	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_rec.no_parents.2.db &> obs
check obs exp
rm obs exp

###################################################################
# 9. Test with two families lacking parents, requiring two kindreds.
###################################################################
echo "    auto_rec.t9...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
WARNING: auto-recessive called on family 2 where no affected has parents
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_rec.no_parents.2.db &> obs
check obs exp
rm obs exp

###################################################################
# 10. Test with three family lacking parents, requiring one kindred.
###################################################################
echo "    auto_rec.t10...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
WARNING: auto-recessive called on family 3 where no affected has parents
WARNING: auto-recessive called on family 2 where no affected has parents
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
SYCE1	chr10	135369531	135369532	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/C,T/C,C/C	3_kid	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_rec.no_parents.3.db &> obs
check obs exp
rm obs exp

###################################################################
# 11. Test with three family lacking parents, requiring two kindreds.
###################################################################
echo "    auto_rec.t11...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
WARNING: auto-recessive called on family 3 where no affected has parents
WARNING: auto-recessive called on family 2 where no affected has parents
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/T,C/T,T/T	2_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/C,T/C,C/C	2_kid	2" > exp

gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_rec.no_parents.3.db &> obs
check obs exp
rm obs exp

###################################################################
# 12. Test with three family lacking parents, requiring two kindreds.
###################################################################
echo "    auto_rec.t12...\c"
echo "WARNING: auto-recessive called on family 1 where no affected has parents
WARNING: auto-recessive called on family 3 where no affected has parents
WARNING: auto-recessive called on family 2 where no affected has parents" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    test.auto_rec.no_parents.3.db &> obs
check obs exp
rm obs exp


###################################################################
# 13. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_rec.t13...\c"
echo "WARNING: no affecteds in family 1
WARNING: no affecteds in family 3
WARNING: no affecteds in family 2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_rec.no_parents.4.db &> obs
check obs exp
rm obs exp

###################################################################
# 14. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_rec.t14...\c"
echo "WARNING: no affecteds in family 1
WARNING: no affecteds in family 3
WARNING: no affecteds in family 2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
    test.auto_rec.no_parents.4.db &> obs
check obs exp
rm obs exp

###################################################################
# 15. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_rec.t15...\c"
echo "WARNING: no affecteds in family 1
WARNING: no affecteds in family 3
WARNING: no affecteds in family 2" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
    test.auto_rec.no_parents.4.db &> obs
check obs exp
rm obs exp


###################################################################
# 16. Test with three family lacking parents and no one affected.
###################################################################
echo "    auto_rec.t16...\c"
echo "WARNING: no affecteds in family 3
WARNING: no affecteds in family 2
gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	2	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	C/T,C/T,T/T	1_kid	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/C,T/C,C/C	1_kid	1" > exp
gemini autosomal_recessive  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
    test.auto_rec.no_parents.5.db &> obs
check obs exp
rm obs exp

###################################################################
# 17. Test with --families
###################################################################
echo "    auto_rec.t17...\c"
#NOTE: don't include unaffecteds by default.
echo "WARNING: no affecteds in family 3
WARNING: no affecteds in family 2" > exp
gemini autosomal_recessive  \
    --columns "chrom, start, end, ref, alt, impact" \
    --families 3 test.auto_rec.no_parents.5.db &> obs
check obs exp
rm obs exp
