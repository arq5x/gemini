###########################################################################################
#1. Test loading a file without using bp-res GERP scores
###########################################################################################
gemini load -v test1.snpeff.vcf --skip-gerp-bp \
	        nogerp.db > /dev/null

echo "    gerp.t1...\c"
echo "chr1	30859	30860	G	C	None	None
chr1	69269	69270	A	G	None	5.41355e-97
chr1	69510	69511	A	G	None	5.41355e-97
chr1	874815	874816	C	CT	None	None
chr1	879675	879676	G	A	None	None
chr1	935491	935492	G	T	None	None
chr1	1334051	1334057	CTAGAG	C	None	1.51184e-59" > exp

gemini query -q "select chrom, start, end, ref, alt, gerp_bp_score, gerp_element_pval from variants" \
	nogerp.db > obs
check obs exp
rm obs exp


###########################################################################################
#2. Test loading a file using bp-res GERP scores
###########################################################################################
gemini load -v test1.snpeff.vcf gerp.db > /dev/null

echo "    gerp.t2...\c"
echo "chr1	30859	30860	G	C	0.0	None
chr1	69269	69270	A	G	-0.726000010967	5.41355e-97
chr1	69510	69511	A	G	1.14999997616	5.41355e-97
chr1	874815	874816	C	CT	0.0267999991775	None
chr1	879675	879676	G	A	-3.1099998951	None
chr1	935491	935492	G	T	0.533999979496	None
chr1	1334051	1334057	CTAGAG	C	4.82000017166	1.51184e-59" > exp

gemini query -q "select chrom, start, end, ref, alt, gerp_bp_score, gerp_element_pval from variants" \
	gerp.db > obs
check obs exp
rm obs exp
