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
# 2. Test comp-het with 
###############################################################################
echo "    comp_het.t1...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor_variant	3	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	TTCT|T,TTCT/T,TTCT/TTCT	child_4	1	1_3_7	3
chr1	17729	17730	WASH7P	C	A	splice_acceptor_variant	7	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	C/A,C/A,C/A	child_4	1	1_3_7	3" > exp.$$

gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--allow-unaffected \
	--max-priority 3 \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs.$$
check obs.$$ exp.$$
rm obs.$$ exp.$$

##### change max priority to default.
echo "    comp_het.t2...\c"

touch exp
gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--allow-unaffected \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp
###############################################################################
# 3. Test comp-het
###############################################################################
echo "    comp_het.t3...\c"

echo "chrom	start	end	gene	ref	alt	impact	exac_num_hom_alt	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor_variant	-1	3	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	TTCT|T,TTCT/T,TTCT/TTCT	child_4	1	1_3_7	3
chr1	17729	17730	WASH7P	C	A	splice_acceptor_variant	-1	7	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	C/A,C/A,C/A	child_4	1	1_3_7	3" > exp

gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact,exac_num_hom_alt" \
    --filter "impact_severity = 'HIGH'" \
	--max-priority 3 \
     test.comp_het.db > obs
check obs exp
rm obs exp


echo "    comp_het.t3b...\c"
touch exp

gemini comp_hets \
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
echo "chrom	start	end	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	17362	17366	1	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	TTCT|T,TTCT|T,TTCT|TTCT	child_3	1	1_1_2	1
chr1	17729	17730	2	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	A|C,C|A,C|A	child_3	1	1_1_2	1" > exp
gemini comp_hets \
    --columns "chrom, start, end" \
	--allow-unaffected \
    --families 3 \
    test.comp_het_default.2.db > obs
check obs exp
rm obs exp

######################################################
# 5. Test basic comp-het functionality with --families
######################################################
echo "    comp_het.t5...\c"
echo "chrom	start	end	ref	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	17362	17366	TTCT	T	1	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	TTCT|T,TTCT|T,TTCT|TTCT	child_3	1	1_1_2	1
chr1	17729	17730	C	A	2	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	A|C,C|A,C|A	child_3	1	1_1_2	1" > exp
gemini comp_hets \
    --columns "chrom, start, end, ref, alt" \
    --families 3 \
    test.comp_het_default.2.db > obs
check obs exp
rm obs exp

###############################################################################
# 6. Test --min-kindreds
###############################################################################
echo "    comp_het.t6...\c"
touch exp
gemini comp_hets \
	--allow-unaffected \
    --min-kindreds 1 \
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
    --min-kindreds 2 \
	--allow-unaffected \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.4.db > obs
check obs exp
rm obs exp

###############################################################################
# 11. Negative Test with --min-kindreds 2 and --only-affected
###############################################################################
echo "    comp_het.t9...\c"
touch exp
gemini comp_hets \
    --min-kindreds 2 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.6.db > obs
check obs exp
rm obs exp

###############################################################################
# 13. Negative Test with --min-kindreds 3 and --only-affected. 
#     Two different comp_hets in same gene
###############################################################################
echo "    comp_het.t10...\c"
touch exp
gemini comp_hets \
    --min-kindreds 3 \
    --columns "chrom,start,end,gene,ref,alt,impact" \
    test.comp_het_default.7.db > obs
check obs exp
rm obs exp

echo "    comp_het.t11"
echo "chrom	start	end	ref	alt	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	17362	17366	TTCT	T	WASH7P	exon_variant	3	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	TTCT|T,TTCT/T,TTCT/TTCT	child_4	1	1_3_5	1
chr1	17729	17730	C	A	WASH7P	exon_variant	5	4	child_4(child_4;affected;male),dad_4(dad_4;unaffected;male),mom_4(mom_4;unaffected;female)	A|C,C/C,C/A	child_4	1	1_3_5	1" > exp

gemini comp_hets \
	--column "chrom,start,end,ref,alt,gene,impact" \
	test.comp_het.unphase.db > obs
check obs exp
rm obs exp

echo "comp_het.t12"
touch exp
gemini comp_hets --columns "chrom, start, end, ref, alt" from_inheritance.db > obs
check obs exp

echo "comp_het.t13"
echo "chrom	start	end	ref	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	100001	100002	A	G	1	1	dad(dad;unknown;male),mom(mom;unknown;female),akid(akid;unknown;male),bkid(bkid;unknown;male)	A/A,A/G,G|A,A/A		1	1_1_2	1
chr1	100002	100003	A	G	2	1	dad(dad;unknown;male),mom(mom;unknown;female),akid(akid;unknown;male),bkid(bkid;unknown;male)	A/G,A/A,A|G,A/A		1	1_1_2	1" > exp
gemini comp_hets --columns "chrom, start, end, ref, alt" --pattern-only from_inheritance.db > obs
check obs exp


echo "comp_het.t14"
echo "start	end	ref	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
100001	100002	A	G	1	1	dad(dad;unknown;male),mom(mom;unknown;female),akid(akid;unknown;male),bkid(bkid;unknown;male)	A/A,A/G,G|A,A/A		1	1_1_2	1
100002	100003	A	G	2	1	dad(dad;unknown;male),mom(mom;unknown;female),akid(akid;unknown;male),bkid(bkid;unknown;male)	A/G,A/A,A|G,A/A		1	1_1_2	1" > exp
gemini comp_hets --columns "start, end, ref, alt" --pattern-only from_inheritance.db > obs
check obs exp


###############################################################################
# test that --min-gq has an effect.
echo "    comp_het.t15...\c"
rm -f exp
touch exp

gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--allow-unaffected \
	--max-priority 2 \
	--min-gq 80 \
    --filter "impact_severity = 'HIGH'" \
     test.comp_het.db > obs
check obs exp
rm obs exp

# singleton
echo "comp_het.t16"
touch exp
gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--max-priority 1 \
	test.comp_het.singleton.db > obs

check obs exp
rm obs exp

echo "comp_het.t17"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	17362	17366	WASH7P	TTCT	T	splice_acceptor_variant	1	1	child_1(child_1;affected;male)	TTCT/T	child_1	1	1_1_2	2
chr1	17729	17730	WASH7P	C	A	splice_acceptor_variant	2	1	child_1(child_1;affected;male)	C/A	child_1	1	1_1_2	2" > exp
gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--max-priority 2 \
	test.comp_het.singleton.db > obs
check obs exp
rm obs exp


# CHANGE: no longer include dad_4 in output since he is unaffected.
echo "    comp_het.t18...\c"
echo "chrom	start	end	gene	ref	alt	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id	priority
chr1	16976	16977	DDX11L1	G	A	downstream_gene_variant	1	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	G/A,G/A,G/A	child_3	1	2_1_4	2
chr1	17562	17563	DDX11L1	G	A	downstream_gene_variant	4	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	G/A,G/G,G/G	child_3	1	2_1_4	2
chr1	17562	17563	DDX11L1	G	A	downstream_gene_variant	4	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	G/A,G/G,G/G	child_3	1	3_4_8	2
chr1	17745	17746	DDX11L1	A	G	downstream_gene_variant	8	3	child_3(child_3;affected;male),dad_3(dad_3;unaffected;male),mom_3(mom_3;unaffected;female)	A/G,A/G,A/G	child_3	1	3_4_8	2" > exp
gemini comp_hets \
    --column "chrom,start,end,gene,ref,alt,impact" \
	--max-priority 2 \
    --families 3 \
    --gene-where 'gene != ""' \
     test.comp_het.db > obs
check obs exp
rm obs exp
