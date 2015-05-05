DB=test.esp.db

gemini query -q "select chrom, end, ref, alt, in_esp, aaf_esp_ea, aaf_esp_aa, aaf_esp_all, exome_chip from variants" test.esp.db > obs

echo "chr1	17588	C	G	0	None	None	None	0
chr1	213068542	C	CTGAT	1	0.20947419433	0.17327080891	0.197140346673	0
chr1	237663970	TA	T	1	0	0.0363534675615	0.0113715885234	0
chr4	110884488	G	GA	1	0	0.0121951219512	0.00415401821377	0" > exp

echo "esp.zero.t1"
check obs exp 
