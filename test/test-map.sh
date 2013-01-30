################################################################################
#1. Test gms mappability
################################################################################
echo "    map.t1...\c"
echo "None	None	None
0.0	16.1	0.0
None	None	None
None	None	None
None	None	None
None	None	None
None	None	None
None	None	None" > exp
gemini query -q "select gms_illumina, gms_iontorrent, gms_solid from variants" test3.snpeff.db \
 > obs
check obs exp
rm obs exp

################################################################################
#2.Test grc mapability
################################################################################

echo "    map.t2...\c"
echo "None
None
None
None
None
None
None
grc_fix" > exp
gemini query -q "select grc from variants" test3.snpeff.db \
> obs
check obs exp
rm obs exp

################################################################################
