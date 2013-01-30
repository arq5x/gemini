################################################################################
#1. Test genome annotations: rmsk
################################################################################
echo "    genome.t1...\c"
echo "Simple_repeat_Simple_repeat_(TC)n;trf;LTR_ERVL-MaLR_MLT1A
None
None
trf
None
None
None
None" > exp
gemini query -q "select rmsk from variants" test3.snpeff.db \
 > obs
check obs exp
rm obs exp

################################################################################
#2.Test genome annotations: in_cpg_island
################################################################################

echo "    genome.t2...\c"
echo "0
0
0
0
0
1
0
0" > exp
gemini query -q "select in_cpg_island from variants" test3.snpeff.db \
> obs
check obs exp
rm obs exp

################################################################################
#3. Test genome annotations: in_segdup
################################################################################
echo "    genome.t3...\c"
echo "1
1
1
0
0
0
0
0" > exp
gemini query -q "select in_segdup from variants" test3.snpeff.db \
 > obs
check obs exp
rm obs exp

################################################################################
#4.Test genome annotations: is_conserved
################################################################################

echo "    genome.t4...\c"
echo "0
0
0
0
0
0
1
1" > exp
gemini query -q "select is_conserved from variants" test3.snpeff.db \
> obs
check obs exp
rm obs exp

################################################################################
#5.Test genome annotations: cyto_band
################################################################################

echo "    genome.t5...\c"
echo "chr1p36.33
chr1p36.33
chr1p36.33
chr1p36.33
chr1p36.33
chr1p36.33
chr1p36.33
chr1p35.2" > exp
gemini query -q "select cyto_band from variants" test3.snpeff.db \
> obs
check obs exp
rm obs exp
################################################################################
#6.Test genome annotations: recomb rate
################################################################################

echo "    genome.t6...\c"
echo "2.981822
2.082414
2.082414
3.162755
3.162755
0.952858
1.84494
0.225399" > exp
gemini query -q "select recomb_rate from variants" test3.snpeff.db \
> obs
check obs exp
rm obs exp

################################################################################