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
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
WDR37	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	1
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	1
ASAH2C	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	C/C,C/C,C/T	3_kid	1
SPRN	G	A	intron	LOW	4	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	G/G,G/G,G/A	1_kid	2
SPRN	G	A	intron	LOW	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	G/G,G/G,G/A	2_kid	2
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    test.de_novo.db > obs
check obs exp
rm obs exp


###################################################################
# 2. Test de_novo with filter
###################################################################
echo "    de_novo.t2...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
WDR37	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	1" > exp
gemini de_novo  \
	--allow-unaffected \
    --columns "gene, ref, alt, impact, impact_severity" \
    --filter "impact_severity = 'HIGH'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test de_novo with filter and minimum depth requirement
###################################################################
echo "    de_novo.t3...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	1" > exp
gemini de_novo  \
	--allow-unaffected \
    --columns "gene, ref, alt, impact, impact_severity" \
    --filter "impact_severity = 'MED'" \
    -d 40 \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test de_novo with filter with min-kindreds
###################################################################
echo "    de_novo.t4...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	2
ASAH2C	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	C/C,C/C,C/T	3_kid	2
SPRN	G	A	intron	LOW	4	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	G/G,G/G,G/A	1_kid	2
SPRN	G	A	intron	LOW	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	G/G,G/G,G/A	2_kid	2
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
	--allow-unaffected \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 5. Test de_novo with filter with min-kindreds 2 and filter
###################################################################
echo "    de_novo.t5...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	2
ASAH2C	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	C/C,C/C,C/T	3_kid	2
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    --min-kindreds 2 \
	--allow-unaffected \
    --filter "impact_severity != 'LOW'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test de_novo with filter with min-kindreds 3 and filter
###################################################################
echo "    de_novo.t6...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    --min-kindreds 3 \
	--allow-unaffected \
    --filter "impact_severity != 'LOW'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 7. Test de_novo with filter with min-kindreds 1 and filter
###################################################################
echo "    de_novo.t7...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	2
ASAH2C	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	C/C,C/C,C/T	3_kid	2
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3
WDR37	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid(1_kid;affected)	T/T,T/T,T/C	1_kid	1" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    --min-kindreds 1 \
	--allow-unaffected \
    --filter "impact_severity != 'LOW'" \
    test.de_novo.db > obs
check obs exp
rm obs exp

###################################################################
# 8. Test de_novo without --only-affected
###################################################################
echo "    de_novo.t8...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
WDR37	T	C	stop_loss	HIGH	1	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid1(1_kid1;unaffected),1_kid2(1_kid2;affected)	T/T,T/T,T/C,T/C	1_kid2	1
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	1
ASAH2C	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	C/C,C/C,C/T	3_kid	1
SPRN	G	A	intron	LOW	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	G/G,G/G,G/A	2_kid	1
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid1(1_kid1;unaffected),1_kid2(1_kid2;affected)	T/T,T/T,T/T,T/C	1_kid2	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
	--allow-unaffected \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp

###################################################################
# 9. Test de_novo with --only-affected
###################################################################
echo "    de_novo.t9...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	1
ASAH2C	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	C/C,C/C,C/T	3_kid	1
SPRN	G	A	intron	LOW	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	G/G,G/G,G/A	2_kid	1
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid1(1_kid1;unaffected),1_kid2(1_kid2;affected)	T/T,T/T,T/T,T/C	1_kid2	3
SYCE1	T	C	non_syn_coding	MED	5	3	3_dad(3_dad;unaffected),3_mom(3_mom;unaffected),3_kid(3_kid;affected)	T/T,T/T,T/C	3_kid	3
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	3" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp


###################################################################
# 10. Test de_novo with --only-affected and --families
###################################################################
echo "    de_novo.t10...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	1
SPRN	G	A	intron	LOW	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	G/G,G/G,G/A	2_kid	1
SYCE1	T	C	non_syn_coding	MED	5	1	1_dad(1_dad;unaffected),1_mom(1_mom;unaffected),1_kid1(1_kid1;unaffected),1_kid2(1_kid2;affected)	T/T,T/T,T/T,T/C	1_kid2	2
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	2" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    --families 1,2 \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp

###################################################################
# 11. Test de_novo with --only-affected and --families
###################################################################
echo "    de_novo.t11...\c"
echo "gene	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	C	T	non_syn_coding	MED	2	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	C/C,C/C,C/T	2_kid	1
SPRN	G	A	intron	LOW	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	G/G,G/G,G/A	2_kid	1
SYCE1	T	C	non_syn_coding	MED	5	2	2_dad(2_dad;unaffected),2_mom(2_mom;unaffected),2_kid(2_kid;affected)	T/T,T/T,T/C	2_kid	1" > exp
gemini de_novo  \
    --columns "gene, ref, alt, impact, impact_severity" \
    --families 2 \
    test.de_novo.affected.and.unaffected.db > obs
check obs exp
rm obs exp
