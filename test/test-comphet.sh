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
    test.comp_het_default.db > obs
echo "chrom	start	end	ref	alt	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	TTCT	T	WASH7P	splice_acceptor	1	0	child_1(unknown),child_2(unknown),dad_2(unknown),mom_2(unknown),dad_1(unknown),mom_1(unknown),child_3(unknown),dad_3(unknown),mom_3(unknown),child_4(unknown),dad_4(unknown),mom_4(unknown)	TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|T,TTCT|T,TTCT|TTCT	child_4	1	1
chr1	17729	17730	C	A	WASH7P	splice_acceptor	2	0	child_1(unknown),child_2(unknown),dad_2(unknown),mom_2(unknown),dad_1(unknown),mom_1(unknown),child_3(unknown),dad_3(unknown),mom_3(unknown),child_4(unknown),dad_4(unknown),mom_4(unknown)	C|C,C|C,C|C,C|C,C|A,C|C,C|A,C|C,C|C,A|C,C|A,C|A	child_4	1	1" > exp
check obs exp
rm obs exp

###############################################################################
# 2. Test comp-het with --ignore-phasing
###############################################################################
echo "    comp_het.t2...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	3	4	dad_4(dad;unaffected),mom_4(mom;unaffected),child_4(child;affected)	TTCT/T,TTCT/TTCT,TTCT/T	child_4	2	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	7	4	dad_4(dad;unaffected),mom_4(mom;unaffected),child_4(child;affected)	C/A,C/A,C/A	child_4	2	1
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	3	4	dad_4(dad;unaffected),mom_4(mom;unaffected),child_4(child;affected)	TTCT/T,TTCT/TTCT,TTCT/T	dad_4	2	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	7	4	dad_4(dad;unaffected),mom_4(mom;unaffected),child_4(child;affected)	C/A,C/A,C/A	dad_4	2	1" > exp
gemini comp_hets --ignore-phasing \
    --column "chrom,start,end,gene,ref,alt,impact" \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp

###############################################################################
# 3. Test comp-het for --only-affected
###############################################################################
echo "    comp_het.t3...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id" > exp
gemini comp_hets --ignore-phasing \
    --only-affected \
    --column "chrom,start,end,gene,ref,alt,impact" \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp

###############################################################################
# 4. Test basic comp-het functionality with --families
###############################################################################
echo "    comp_het.t4...\c"
echo "chrom	start	end	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	1	3	dad_3(dad;unaffected),mom_3(mom;unaffected),child_3(child;affected)	TTCT|T,TTCT|TTCT,TTCT|T	child_3	2	1
chr1	17729	17730	2	3	dad_3(dad;unaffected),mom_3(mom;unaffected),child_3(child;affected)	C|A,C|A,A|C	child_3	2	1
chr1	17362	17366	1	3	dad_3(dad;unaffected),mom_3(mom;unaffected),child_3(child;affected)	TTCT|T,TTCT|TTCT,TTCT|T	dad_3	2	1
chr1	17729	17730	2	3	dad_3(dad;unaffected),mom_3(mom;unaffected),child_3(child;affected)	C|A,C|A,A|C	dad_3	2	1" > exp
gemini comp_hets \
    --ignore-phasing \
    --columns "chrom, start, end" \
    --families 3 \
    test.comp_het_default.2.db > obs
check obs exp
rm obs exp

###############################################################################
# 5. Test basic comp-het functionality with --families and --only-affected
###############################################################################
echo "    comp_het.t5...\c"
echo "chrom	start	end	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id" > exp
gemini comp_hets \
    --ignore-phasing \
    --columns "chrom, start, end" \
    --families 3 \
    --only-affected \
    test.comp_het_default.2.db > obs
check obs exp
rm obs exp

###############################################################################
# 6. Test --min-kindreds
###############################################################################
echo "    comp_het.t6...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_1	1	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	C|C,C|C,C|A	child_1	1	1" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 1 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 7. Test without --min-kindreds, but default should be 1
###############################################################################
echo "    comp_het.t7...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_1	1	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	C|C,C|C,C|A	child_1	1	1" > exp
gemini comp_hets \
    --ignore-phasing \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 8. Negative Test with --min-kindreds 2
###############################################################################
echo "    comp_het.t8...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 9. Positive Test with --min-kindreds 2
###############################################################################
echo "    comp_het.t9...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_1	2	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	C|C,C|C,C|A	child_1	2	1
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	2	dad_2(dad;unaffected),mom_2(mom;unaffected),child_2(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_2	2	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	2	dad_2(dad;unaffected),mom_2(mom;unaffected),child_2(child;affected)	C|C,C|C,C|A	child_2	2	1" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.5.db > obs
check obs exp
rm obs exp

###############################################################################
# 9. Positive Test with --min-kindreds 2 and --only-affected
###############################################################################
echo "    comp_het.t10...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_1	2	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	C|C,C|C,C|A	child_1	2	1
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	2	dad_2(dad;unaffected),mom_2(mom;unaffected),child_2(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_2	2	1
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	2	dad_2(dad;unaffected),mom_2(mom;unaffected),child_2(child;affected)	C|C,C|C,C|A	child_2	2	1" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --only-affected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.5.db > obs
check obs exp
rm obs exp

###############################################################################
# 11. Negative Test with --min-kindreds 2 and --only-affected
###############################################################################
echo "    comp_het.t11...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --only-affected \
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
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor	1	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	TTCT|TTCT,TTCT|TTCT,TTCT|T	child_1	1	2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	1	dad_1(dad;unaffected),mom_1(mom;unaffected),child_1(child;affected)	C|C,C|C,C|A	child_1	2	2
chr1	17729	17730	WASH7P	C	A	splice_acceptor	2	2	dad_2(dad;unaffected),mom_2(mom;unaffected),child_2(child;affected)	C|C,C|C,C|A	child_2	2	1
chr1	27729	27730	WASH7P	C	A	splice_acceptor	4	2	dad_2(dad;unaffected),mom_2(mom;unaffected),child_2(child;affected)	C|C,C|C,C|A	child_2	1	1" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --only-affected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.7.db > obs
check obs exp
rm obs exp

###############################################################################
# 13. Negative Test with --min-kindreds 3 and --only-affected. 
#     Two different comp_hets in same gene
###############################################################################
echo "    comp_het.t13...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 3 \
    --only-affected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.7.db > obs
check obs exp
rm obs exp
