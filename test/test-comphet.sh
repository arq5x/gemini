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
echo "family	sample	comp_het_id	chrom	start	end	ref	alt	gene	impact
0	child_4	1	chr1	17362	17366	TTCT	T	WASH7P	splice_acceptor
0	child_4	1	chr1	17729	17730	C	A	WASH7P	splice_acceptor" > exp
gemini comp_hets \
    --column "chrom,start,end,ref,alt,gene,impact" \
    test.comp_het_default.db > obs
check obs exp
rm obs exp

###############################################################################
# 2. Test comp-het with --ignore-phasing
###############################################################################
echo "    comp_het.t2...\c"
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact
4	child_4	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
4	child_4	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor
4	dad_4	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
4	dad_4	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	alt
3	child_3	1	chr1	17362	17366	WASH7P	T
3	child_3	1	chr1	17729	17730	WASH7P	A
3	dad_3	1	chr1	17362	17366	WASH7P	T
3	dad_3	1	chr1	17729	17730	WASH7P	A" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	alt" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact
1	child_1	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
1	child_1	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact
1	child_1	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
1	child_1	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact
1	child_1	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
1	child_1	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor
2	child_2	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
2	child_2	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact
1	child_1	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
1	child_1	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor
2	child_2	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
2	child_2	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact" > exp
gemini comp_hets \
    --ignore-phasing \
    --min-kindreds 2 \
    --only-affected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.6.db > obs
check obs exp
rm obs exp
