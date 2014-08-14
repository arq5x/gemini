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
             --gt-filter "(gt_types).(phenotype==1).(==HOM_REF).(all)" extended_ped.db > obs
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
             --gt-filter "(gt_types).(phenotype==1).(==HOM_REF).(all)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 6. Test a wildcard filter with a vanilla filter
########################################################################
echo "    wildcard.t6...\c"
echo "chrom	start	end	ref	alt	gts.M10478	gts.M10500
chr10	135336655	135336656	G	A	A/A	./." > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==2) from variants" \
             --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) and gts.M10478 =='A/A'" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 7. Test specific wildcard with genotype filter using ANY
########################################################################
echo "    wildcard.t7...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M128215
chr10	48003991	48003992	C	T	T/T	C/C
chr10	126678091	126678092	G	A	G/G	G/A
chr10	135210790	135210791	T	C	T/T	T/T
chr10	135369531	135369532	T	C	T/T	T/T
chr16	72057434	72057435	C	T	C/T	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==1) from variants" \
             --gt-filter "(gt_types).(phenotype==1).(==HOM_REF).(any)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 8. Test specific wildcard with genotype filter using NONE
########################################################################
echo "    wildcard.t8...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M128215
chr10	1142207	1142208	T	C	C/C	C/C
chr10	52004314	52004315	T	C	./.	C/C
chr10	52497528	52497529	G	C	./.	./.
chr10	135336655	135336656	G	A	./.	A/A" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==1) from variants" \
             --gt-filter "(gt_types).(phenotype==1).(==HOM_REF).(none)" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 9. Test specific wildcard with genotype filter using ANY and *
########################################################################
echo "    wildcard.t9...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	48003991	48003992	C	T	T/T	C/T	C/T	C/C
chr10	126678091	126678092	G	A	G/G	G/G	G/G	G/A
chr10	135369531	135369532	T	C	T/T	T/C	T/C	T/T
chr16	72057434	72057435	C	T	C/T	C/C	C/C	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HET).(any)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 10. Test specific wildcard with genotype filter using COUNT > 0 (should be same as any)
########################################################################
echo "    wildcard.t10...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	48003991	48003992	C	T	T/T	C/T	C/T	C/C
chr10	126678091	126678092	G	A	G/G	G/G	G/G	G/A
chr10	135369531	135369532	T	C	T/T	T/C	T/C	T/T
chr16	72057434	72057435	C	T	C/T	C/C	C/C	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HET).(count>0)" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 11. Test specific wildcard with genotype filter using COUNT >= 2 
########################################################################
echo "    wildcard.t11...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	48003991	48003992	C	T	T/T	C/T	C/T	C/C
chr10	135369531	135369532	T	C	T/T	T/C	T/C	T/T" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HET).(count>=2)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 12. Test specific wildcard with genotype filter using COUNT == 0  
########################################################################
echo "    wildcard.t12...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	1142207	1142208	T	C	C/C	C/C	C/C	C/C
chr10	52004314	52004315	T	C	./.	./.	C/C	C/C
chr10	52497528	52497529	G	C	./.	C/C	C/C	./.
chr10	135210790	135210791	T	C	T/T	C/C	C/C	T/T
chr10	135336655	135336656	G	A	./.	A/A	./.	A/A" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HET).(count==0)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 13. Test specific wildcard with genotype filter using COUNT != 2  
########################################################################
echo "    wildcard.t13...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	1142207	1142208	T	C	C/C	C/C	C/C	C/C
chr10	48003991	48003992	C	T	T/T	C/T	C/T	C/C
chr10	126678091	126678092	G	A	G/G	G/G	G/G	G/A
chr10	135369531	135369532	T	C	T/T	T/C	T/C	T/T
chr16	72057434	72057435	C	T	C/T	C/C	C/C	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HOM_ALT).(count!=2)" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 14. Test specific wildcard with genotype filter using COUNT != 2 and specific filter  
########################################################################
echo "    wildcard.t14...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	126678091	126678092	G	A	G/G	G/G	G/G	G/A" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HOM_ALT).(count!=2) and gts.M10478 =='G/G'" extended_ped.db > obs
check obs exp
rm obs exp


########################################################################
# 15. Test specific wildcard with genotype filter using different spacing
########################################################################
echo "    wildcard.t15...\c"
echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	126678091	126678092	G	A	G/G	G/G	G/G	G/A" > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==    HOM_ALT).(  count   !=2) and gts.M10478 =='G/G'" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 16. Test syntax failure  
########################################################################
echo "    wildcard.t16...\c"
echo "Unsupported wildcard operation: (). Exiting." > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HOM_ALT).()" extended_ped.db 2> obs
check obs exp
rm obs exp

########################################################################
# 17. Test syntax failure  
########################################################################
echo "    wildcard.t17...\c"
echo "Unsupported wildcard operation: (amy). Exiting." > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HOM_ALT).(amy)" extended_ped.db 2> obs
check obs exp
rm obs exp

########################################################################
# 18. Test syntax failure  
########################################################################
echo "    wildcard.t18...\c"
echo "Wildcard filter should consist of 4 elements. Exiting." > exp
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
             --gt-filter "(gt_types).(*).(==HOM_ALT)" extended_ped.db 2> obs
check obs exp
rm obs exp

########################################################################
# 19. Test multiple wildcard on same column
########################################################################
echo "    wildcard.t19...\c"
echo "chrom	start	end	ref	alt	gene	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	1142207	1142208	T	C	WDR37	C/C	C/C	C/C	C/C
chr10	48003991	48003992	C	T	ASAH2C	T/T	C/T	C/T	C/C
chr10	52004314	52004315	T	C	ASAH2	./.	./.	C/C	C/C
chr10	52497528	52497529	G	C	ASAH2B	./.	C/C	C/C	./.
chr10	135336655	135336656	G	A	SPRN	./.	A/A	./.	A/A" > exp
gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(*) from variants" \
             --gt-filter "((gt_types).(phenotype==1).(!=HOM_REF).(count>=1) and (gt_types).(phenotype==2).(!=HOM_REF).(count>=1))" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 20. Test multiple wildcard on same column
########################################################################
echo "    wildcard.t20...\c"
echo "chrom	start	end	ref	alt	gene	gts.M10475	gts.M10478	gts.M10500	gts.M128215" > exp
gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(*) from variants" \
             --gt-filter "((gt_types).(phenotype==1).(!=HOM_REF).(all) and (gt_types).(phenotype==2).(==HOM_REF).(all))" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 21. Test multiple wildcard on same column
########################################################################
echo "    wildcard.t21...\c"
echo "chrom	start	end	ref	alt	gene	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	126678091	126678092	G	A	CTBP2	G/G	G/G	G/G	G/A
chr16	72057434	72057435	C	T	DHODH	C/T	C/C	C/C	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(*) from variants" \
             --gt-filter "((gt_types).(phenotype==1).(!=HOM_REF).(count>=1) and (gt_types).(phenotype==2).(==HOM_REF).(all))" extended_ped.db > obs
check obs exp
rm obs exp

########################################################################
# 22. Test AND 
########################################################################
echo "    wildcard.t22...\c"
echo "chrom	start	end	ref	alt	gene	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	126678091	126678092	G	A	CTBP2	G/G	G/G	G/G	G/A
chr16	72057434	72057435	C	T	DHODH	C/T	C/C	C/C	C/C" > exp
gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(*) from variants" \
             --gt-filter "((gt_types).(phenotype==1).(!=HOM_REF).(count>=1) AND (gt_types).(phenotype==2).(==HOM_REF).(all))" extended_ped.db > obs
check obs exp
rm obs exp
