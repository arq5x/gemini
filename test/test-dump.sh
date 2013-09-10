####################################################################
# 1. Test the TFAM dump
####################################################################
echo "    dump.t01...\c"
echo "None M10475 None None None None
None M10478 None None None None
None M10500 None None None None
None M128215 None None None None" > exp
gemini dump --tfam test4.snpeff.db > obs
check obs exp
rm obs exp
