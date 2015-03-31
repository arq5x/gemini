source ./check.sh
################################################################################
#1. Test ESP annotations
################################################################################
echo "0	None	None	None	0
1	0.0457067757009	0.003663003663	0.0306466729147	0
1	0.88742933156	0.544101123596	0.759766033006	0
1	0.136133389616	0.0213645761544	0.0925392670157	0
0	None	None	None	0
1	0.00511985105888	0.000227686703097	0.00346527029108	1
0	None	None	None	0
0	None	None	None	0
0	None	None	None	0" > exp
gemini query -q "select in_esp, aaf_esp_ea, aaf_esp_aa, aaf_esp_all, exome_chip from variants" test2.snpeff.db \
 > obs
check obs exp "pop_metrics.t1"
rm obs exp

################################################################################
#2.Test 1000g annotations
################################################################################

echo "in_1kg	aaf_1kg_all	aaf_1kg_eas	aaf_1kg_sas	aaf_1kg_amr	aaf_1kg_afr	aaf_1kg_eur
0	None	None	None	None	None	None
1	0.0189696	0.003	0.0153	0.036	0.0015	0.0497
0	None	None	None	None	None	None
1	0.0367412	0.001	0.045	0.0418	0.0045	0.1034
1	0.688099	0.876	0.8098	0.6254	0.407	0.7942
1	0.00279553	0	0.0041	0.0072	0	0.005
1	0.000599042	0	0.0031	0	0	0
0	None	None	None	None	None	None
1	1	1	1	1	1	1" > exp
gemini query --header -q "select in_1kg, aaf_1kg_all, aaf_1kg_eas, aaf_1kg_sas, aaf_1kg_amr, aaf_1kg_afr, aaf_1kg_eur from variants" test2.snpeff.db \
> obs
check obs exp "pop_metrics.t2"
rm obs exp

################################################################################
#3. Test dbsnp and rsids
################################################################################

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
check obs exp "pop_metrics.t3"
rm obs exp
