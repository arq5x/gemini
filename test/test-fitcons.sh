#############################################
# 1. Test ExAC allele frequencies
#############################################
source ./check.sh

gemini query -q "select fitcons from variants" test.exac.db > obs

echo "0.706548
0.088506
0.646311
0.455138
0.053691
0.615465
0.497415
0.706548" > exp

check obs exp fitcons.t01
rm obs exp


#############################################

gemini query -q "select fitcons from variants" test.amend.db  > obs
echo "0.156188
0.553676
0.061011
0.080785
0.722319
0.069865
0.056701
0.487112
0.732398" > exp

check obs exp fitcons.t02
rm obs exp


#############################################


gemini query -q "select fitcons from variants" test.query.vep.db | awk '(NR % 50 == 0)'  > obs
echo "0.095506
0.099367
0.121938
0.088506
0.166803
0.118714
0.078448
0.156188
0.055982
0.056701
0.056701
0.114469
0.069865
0.283894
None
0.078448
None" > exp

check obs exp fitcons.t03
