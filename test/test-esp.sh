DB=test.esp.db
source ./check.sh

gemini query -q "select in_esp, aaf_esp_ea, aaf_esp_aa, aaf_esp_all, exome_chip from variants" test.esp.db > obs

echo "0	None	None	None	0
0	None	None	None	0
0	None	None	None	0
1	0	0.0714285714286	0.0714285714286	0" > exp

check obs exp "esp.zero.t1"
