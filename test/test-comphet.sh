check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}
export -f check
###############################################################################
# 1. Test basic comp-het functionality (for phased genotypes)
###############################################################################
echo "    comp_het.t1...\c"
# NOTE the families all have id=0 in the db, so they all show up.
gemini comp_hets \
    --column "chrom,start,end,ref,alt,gene,impact" \
	--allow-unaffected \
	--lenient \
    test.comp_het_default.db > obs
echo "chrom	start	end	ref	alt	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	TTCT	T	WASH7P	splice_acceptor	1	0	child_1(child_1;unknown),child_2(child_2;unknown),dad_2(dad_2;unknown),mom_2(mom_2;unknown),dad_1(dad_1;unknown),mom_1(mom_1;unknown),child_3(child_3;unknown),dad_3(dad_3;unknown),mom_3(mom_3;unknown),child_4(child_4;unknown),dad_4(dad_4;unknown),mom_4(mom_4;unknown)	TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|T,TTCT|T,TTCT|TTCT		1	1_1_2
chr1	17729	17730	C	A	WASH7P	splice_acceptor	2	0	child_1(child_1;unknown),child_2(child_2;unknown),dad_2(dad_2;unknown),mom_2(mom_2;unknown),dad_1(dad_1;unknown),mom_1(mom_1;unknown),child_3(child_3;unknown),dad_3(dad_3;unknown),mom_3(mom_3;unknown),child_4(child_4;unknown),dad_4(dad_4;unknown),mom_4(mom_4;unknown)	C|C,C|C,C|C,C|C,C|A,C|C,C|A,C|C,C|C,A|C,C|A,C|A		1	1_1_2" > exp

check obs exp
rm obs exp

###############################################################################
# 2. Test comp-het with --ignore-phasing
###############################################################################
# CHANGE: no longer include dad_4 in output since he is unaffected.
echo "    comp_het.t2...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	3	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	TTCT|T,TTCT/T,TTCT/TTCT	child_4	1	1_3_7
chr1	17729	17730	WASH7P	C	A	splice_acceptor	7	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	C/A,C/A,C/A	child_4	1	1_3_7" > exp

gemini comp_hets \
	--ignore-phasing \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--allow-unaffected \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp

###############################################################################
# 3. Test comp-het for --only-affected
###############################################################################
echo "    comp_het.t3...\c"
touch exp
gemini comp_hets --ignore-phasing \
    --column "chrom,start,end,gene,ref,alt,impact" \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp

###############################################################################
# 4. Test basic comp-het functionality with --families
###############################################################################
echo "    comp_het.t4...\c"
#CHANGE: only output child 3 because dad_3 is not affected.
echo "chrom	start	end	gene	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	T	1	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	TTCT|T,TTCT|T,TTCT|TTCT	child_3	1	1_1_2
chr1	17729	17730	WASH7P	A	2	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	A|C,C|A,C|A	child_3	1	1_1_2" > exp
gemini comp_hets \
    --ignore-phasing \
    --columns "chrom, start, end" \
	--allow-unaffected \
    --families 3 \
    test.comp_het_default.2.db > obs
check obs exp
rm obs exp

###############################################################################
# 5. Test basic comp-het functionality with --families and --only-affected
###############################################################################
echo "    comp_het.t5...\c"
touch exp
gemini comp_hets \
    --ignore-phasing \
    --columns "chrom, start, end" \
    --families 3 \
    test.comp_het_default.2.db > obs
check obs exp
rm obs exp

###############################################################################
# 6. Test --min-kindreds
###############################################################################
echo "    comp_het.t6...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_1	1	1_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	C|A,C|C,C|C	child_1	1	1_1_2" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 1 \
	--allow-unaffected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 7. Test without --min-kindreds, but default should be 1
###############################################################################
echo "    comp_het.t7...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_1	1	1_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	C|A,C|C,C|C	child_1	1	1_1_2" > exp
gemini comp_hets \
    --ignore-phasing \
	--allow-unaffected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 8. Negative Test with --min-kindreds 2
###############################################################################
echo "    comp_het.t8...\c"
touch exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
	--allow-unaffected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 9. Positive Test with --min-kindreds 2
###############################################################################
echo "    comp_het.t9...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_1	2	1_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	C|A,C|C,C|C	child_1	2	1_1_2
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	2	child_2(child_2;affected;male),dad_2(dad_2;unaffected;male),mom_2(mom_2;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_2	2	1_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	2	child_2(child_2;affected;male),dad_2(dad_2;unaffected;male),mom_2(mom_2;unaffected;female)	C|A,C|C,C|C	child_2	2	1_1_2" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
	--allow-unaffected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.5.db > obs
check obs exp
rm obs exp

###############################################################################
# 9. Positive Test with --min-kindreds 2 and --only-affected
###############################################################################
echo "    comp_het.t10...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_1	2	1_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	C|A,C|C,C|C	child_1	2	1_1_2
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	2	child_2(child_2;affected;male),dad_2(dad_2;unaffected;male),mom_2(mom_2;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_2	2	1_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	2	child_2(child_2;affected;male),dad_2(dad_2;unaffected;male),mom_2(mom_2;unaffected;female)	C|A,C|C,C|C	child_2	2	1_1_2" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.5.db > obs
check obs exp
rm obs exp

###############################################################################
# 11. Negative Test with --min-kindreds 2 and --only-affected
###############################################################################
echo "    comp_het.t11...\c"
touch exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.6.db > obs
check obs exp
rm obs exp

###############################################################################
# 12. Positive Test with --min-kindreds 2 and --only-affected. 
#     Two different comp_hets in same gene
###############################################################################
echo "    comp_het.t12...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	2	child_2(child_2;affected;male),dad_2(dad_2;unaffected;male),mom_2(mom_2;unaffected;female)	C|A,C|C,C|C	child_2	2	1_2_4
chr1	27729	27730	WASH7P	C	A	splice_acceptor	4	2	child_2(child_2;affected;male),dad_2(dad_2;unaffected;male),mom_2(mom_2;unaffected;female)	C|A,C|C,C|C	child_2	2	1_2_4
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	TTCT|T,TTCT|TTCT,TTCT|TTCT	child_1	2	2_1_2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	child_1(child_1;affected;male),dad_1(dad_1;unaffected;male),mom_1(mom_1;unaffected;female)	C|A,C|C,C|C	child_1	2	2_1_2" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.7.db > obs
check obs exp
rm obs exp

###############################################################################
# 13. Negative Test with --min-kindreds 3 and --only-affected. 
#     Two different comp_hets in same gene
###############################################################################
echo "    comp_het.t13...\c"
touch exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 3 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.7.db > obs
check obs exp
rm obs exp

echo "    comp_het.t14"
echo "chrom	start	end	ref	alt	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	TTCT	T	WASH7P	splice_acceptor	3	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	TTCT|T,TTCT/T,TTCT/TTCT	child_4	1	1_3_5
chr1	17729	17730	C	A	WASH7P	splice_acceptor	5	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	A|C,C/C,C/A	child_4	1	1_3_5" > exp
gemini comp_hets \
	--column "chrom,start,end,ref,alt,gene,impact" \
	test.comp_het.unphase.db > obs
check obs exp
rm obs exp
