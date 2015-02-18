check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}
export -f check
####################################################################
# 1. Test set_somatic for structural variants
####################################################################
echo -e "    fusions.t01...\c"
echo "H_LS-E2-A14P-01A-31D-A19H-09	T/[3:176909982[T	0.235772357724	29	123	H_LS-E2-A14P-10A-01D-A19H-09	T/T	0.0	0	90	chr3	178906029	178906030	T	[3:176909982[T	PIK3CA
H_LS-E2-A14P-01A-31D-A19H-09	G/[3:178906030[G	0.235772357724	29	123	H_LS-E2-A14P-10A-01D-A19H-09	G/G	0.0	0	90	chr3	176909981	176909982	G	[3:178906030[G	TBL1XR1
H_LS-E2-A14P-01A-31D-A19H-09	A/A]X:6146702]	0.253968253968	16	63	H_LS-E2-A14P-10A-01D-A19H-09	A/A	0.0	0	112	chr18	64284795	64284796	A	A]X:6146702]	None
H_LS-E2-A14P-01A-31D-A19H-09	C/C]18:64284796]	0.253968253968	16	63	H_LS-E2-A14P-10A-01D-A19H-09	C/C	0.0	0	112	chrX	6146701	6146702	C	C]18:64284796]	NLGN4X
Identified and set 4 somatic mutations" > exp

gemini set_somatic test.fusions.db > obs

check obs exp
rm obs exp

####################################################################
# 2. Test fusions
####################################################################
echo -e "    fusions.t02...\c"

echo "chr3	176909953	176909982	chr3	178906001	178906030	1233	9.58	-	+	complex	TBL1XR1	PIK3CA	LUMPY	PE	0	H_LS-E2-A14P-01A-31D-A19H-09" > exp

gemini fusions test.fusions.db > obs

check obs exp
rm obs exp
