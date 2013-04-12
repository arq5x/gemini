################################################################################
#1. Test ESP annotations
################################################################################
echo "    pop_metrics.t1...\c"
echo "0	None	None	None	0
1	0.045707	0.003663	0.030647	0
1	0.112571	0.455899	0.240234	0
1	0.136133	0.021365	0.092539	0
0	None	None	None	0
1	0.00512	0.000228	0.003465	1
0	None	None	None	0
0	None	None	None	0
0	None	None	None	0" > exp
gemini query -q "select in_esp, aaf_esp_ea, aaf_esp_aa, aaf_esp_all, exome_chip from variants" test2.snpeff.db \
 > obs
check obs exp
rm obs exp

################################################################################
#2.Test 1000g annotations
################################################################################

echo "    pop_metrics.t2...\c"
echo "0	None	None	None	None	None
0	None	None	None	None	None
1	0.65	0.87	0.65	0.33	0.7
0	None	None	None	None	None
0	None	None	None	None	None
1	0.0032	None	0.01	None	0.01
0	None	None	None	None	None
0	None	None	None	None	None
1	1	0.99	1	1	0.99" > exp
gemini query -q "select in_1kg, aaf_1kg_all, aaf_1kg_asn, aaf_1kg_amr, aaf_1kg_afr, aaf_1kg_eur from variants" test2.snpeff.db \
> obs
check obs exp
rm obs exp

################################################################################
#3. Test dbsnp and rsids
################################################################################

echo "    pop_metrics.t3...\c"
echo "1	rs201219564
1	rs140739101
1	rs75062661
1	rs200505207
1	rs200676709
1	rs41285790
1	rs201898716
0	None
1	rs6672356" > exp
gemini query -q "select in_dbsnp, rs_ids from variants" test2.snpeff.db > obs
check obs exp
rm obs exp
