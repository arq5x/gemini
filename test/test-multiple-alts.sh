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
# 1. Test accessing individual genotype alleles
####################################################################

gemini query --header -q "select (gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*) from variants" --gt-filter "(gt_types).(*).(==2).(all)" test.mad.db > obs2
gemini query --header -q "select (gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*) from variants" --gt-filter "(gt_types).(*).(==1).(all)" test.mad.db > obs1
gemini query --header -q "select (gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*) from variants" --gt-filter "(gt_types).(*).(==0).(all)" test.mad.db > obs0
# this would be empty.
#gemini query --header -q "(gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*) from variants" --gt-filter "(gt_types).(*).(==3).(all)" test.mad.db > obs3

echo "    multiple-alts.t01...\c"
rm -f exp && touch exp
# should be empty because 1 is hets
awk 'NR > 1 { split($1,a, "/"); if($1 == "." || a[1] == a[2]) { print $0 }  }' obs1 > obs
check obs exp

echo "    multiple-alts.t02...\c"
# this should not filter any as 0 is homozygous ref:
awk '(NR == 1) { print $0} NR > 1 { split($1,a, "/"); if($1 == "." || a[1] == a[2]) { print $0 }  }' obs0 > obs
check obs obs0


echo "    multiple-alts.t03...\c"

echo "0" > exp
# unknowns should have 0 for alt depth.
awk 'BEGIN{s=0} $4 > 0 && NR > 1 {s+=1}END{print s}' obs2 > obs
check obs exp

# unknowns should have 0 for and ref depth.
echo "    multiple-alts.t04...\c"
awk 'BEGIN{s=0} $3 > 0 && NR > 1 {s+=1}END{print s}' obs2 > obs
check obs exp

rm -f obs obs1 obs2 obs0 exp
