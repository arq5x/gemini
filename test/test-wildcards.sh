check()
{
  if diff $1 $2; then
      echo ok
  else
      echo fail
  fi
}
export -f check

########################################################################
# 1. Test select * wildcard
########################################################################
echo "    wildcard.t1...\c"
echo "C/C	C/C	C/C	C/C
T/T	C/T	C/T	C/C
./.	./.	C/C	C/C
./.	C/C	C/C	./.
G/G	G/G	G/G	G/A
T/T	C/C	C/C	T/T
./.	A/A	./.	A/A
T/T	T/C	T/C	T/T
C/T	C/C	C/C	C/C" > exp
gemini query -q "select (gts).(*)from variants" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 2. Test select * wildcard with other columns and a header
########################################################################
echo "    wildcard.t2...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	1142207	1142208	T	C	C/C	C/C	C/C	C/C
chr10	48003991	48003992	C	T	T/T	C/T	C/T	C/C
chr10	52004314	52004315	T	C	./.	./.	C/C	C/C
chr10	52497528	52497529	G	C	./.	C/C	C/C	./.
chr10	126678091	126678092	G	A	G/G	G/G	G/G	G/A
chr10	135210790	135210791	T	C	T/T	C/C	C/C	T/T
chr10	135336655	135336656	G	A	./.	A/A	./.	A/A
chr10	135369531	135369532	T	C	T/T	T/C	T/C	T/T
chr16	72057434	72057435	C	T	C/T	C/C	C/C	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*)from variants" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 3. Test specific wildcard with other columns and a header
########################################################################
echo "    wildcard.t3...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M128215
chr10	1142207	1142208	T	C	C/C	C/C
chr10	48003991	48003992	C	T	T/T	C/C
chr10	52004314	52004315	T	C	./.	C/C
chr10	52497528	52497529	G	C	./.	./.
chr10	126678091	126678092	G	A	G/G	G/A
chr10	135210790	135210791	T	C	T/T	T/T
chr10	135336655	135336656	G	A	./.	A/A
chr10	135369531	135369532	T	C	T/T	T/T
chr16	72057434	72057435	C	T	C/T	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==1) from variants" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 4. Test specific wildcard with genotype filter
########################################################################
echo "    wildcard.t4...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M128215
chr10	135210790	135210791	T	C	T/T	T/T
chr10	135369531	135369532	T	C	T/T	T/T" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==1) from variants" \
             --gt-filter "(gt_types).(phenotype==1).(==HOM_REF)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 5. Test multiple select wildcards with genotype filter
########################################################################
echo "    wildcard.t5...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M128215	gt_types.M10475	gt_types.M128215
chr10	135210790	135210791	T	C	T/T	T/T	0	0
chr10	135369531	135369532	T	C	T/T	T/T	0	0" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==1), (gt_types).(phenotype==1) from variants" \
             --gt-filter "(gt_types).(phenotype==1).(==HOM_REF)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 6. Test a wildcard filter with a vanilla filter
########################################################################
echo "    wildcard.t6...\c"
echo "chrom	start	end	ref	alt	gts.M10478	gts.M10500
chr10	135336655	135336656	G	A	A/A	./." > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==2) from variants" \
             --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF) and gts.M10478 =='A/A'" extended_ped.db > obs
check obs exp
rm obs exp





 