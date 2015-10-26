echo "test.vep-extra.t1"
gemini query -q "select impact,vep_impact from variants" test.vep.extra.db > obs
echo "missense_variant	MODERATE" > exp
check obs exp

