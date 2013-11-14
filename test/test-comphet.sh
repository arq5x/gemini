###############################################################################
# 1. Test basic comp-het functionality (for phased genotypes)
###############################################################################
echo "    comp_het.t1...\c"
echo "family	sample	comp_het_id	chrom	start	end	ref	alt	gene	impact
None	child_4	1	chr1	17362	17366	TTCT	T	WASH7P	splice_acceptor
None	child_4	1	chr1	17729	17730	C	A	WASH7P	splice_acceptor" > exp
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
4	dad_4	2	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
4	dad_4	2	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
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
echo "family	sample	comp_het_id	chrom	start	end	gene	ref	alt	impact
4	child_4	1	chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor
4	child_4	1	chr1	17729	17730	WASH7P	C	A	splice_acceptor" > exp
gemini comp_hets --ignore-phasing \
    --only-affected \
    --column "chrom,start,end,gene,ref,alt,impact" \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp






