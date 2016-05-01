echo "test.vep-extra.t1"
gemini query -q "select impact,vep_impact from variants" test.vep.extra.db > obsz
echo "missense_variant	MODERATE
None	None" > expz
check obsz expz

