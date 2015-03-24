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
echo "    genotypes.t01...\c"
echo "./.	./.	./.	./.
G/G	G/G	G/G	G/G
CCT/CCT	CCT/CCT	CCT/C	CCT/CCT
T/C	T/C	T/T	T/T
./.	./.	./.	./.
./.	./.	G/G	G/G
T/T	T/T	T/T	T/T
./.	./.	A/G	A/G
A/A	A/T	A/A	A/A
./.	G/G	G/G	G/G" > exp
gemini query -q "select gts.1094PC0005, gts.1094PC0009, \
				gts.1094PC0012, gts.1094PC0013 \
				from variants" test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 2. Test accessing individual genotype types
####################################################################
echo "    genotypes.t02...\c"
echo "2	2	2	2
0	0	0	0
0	0	1	0
1	1	0	0
2	2	2	2
2	2	3	3
0	0	0	0
2	2	1	1
0	1	0	0
2	0	0	0" > exp
gemini query -q "select gt_types.1094PC0005, gt_types.1094PC0009, \
	                    gt_types.1094PC0012, gt_types.1094PC0013 \
	             from variants" test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 3. Test filtering individual genotype types
####################################################################
echo "    genotypes.t03...\c"
echo "0	0	1	0
2	2	1	1" > exp
gemini query -q "select gt_types.1094PC0005, gt_types.1094PC0009, \
	                    gt_types.1094PC0012, gt_types.1094PC0013 \
	             from variants" \
			 --gt-filter "gt_types.1094PC0012 == HET" \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 4. Test a more complex genotype filter
####################################################################
echo "    genotypes.t04...\c"
echo "chrom	end	ref	alt	gt_types.1094PC0005	gt_types.1094PC0009	gt_types.1094PC0012	gt_types.1094PC0013
chr1	30869	CCT	C	0	0	1	0
chr1	30895	T	C	1	1	0	0
chr1	69511	A	G	2	2	1	1" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_types.1094PC0005, gt_types.1094PC0009, \
	                    gt_types.1094PC0012, gt_types.1094PC0013 \
	             from variants" \
			 --gt-filter "(gt_types.1094PC0012 == HET or \
						   gt_types.1094PC0005 == HET)" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 5.  Test accessing individual genotype depths
####################################################################
echo "    genotypes.t05...\c"
echo "chrom	end	ref	alt	gt_depths.1094PC0005	gt_depths.1094PC0009	gt_depths.1094PC0012	gt_depths.1094PC0013
chr1	30548	T	G	-1	-1	-1	-1
chr1	30860	G	C	7	2	6	4
chr1	30869	CCT	C	8	3	6	5
chr1	30895	T	C	8	3	6	5
chr1	30923	G	T	-1	-1	-1	-1
chr1	69270	A	G	-1	-1	3	6
chr1	69428	T	G	2	79	87	107
chr1	69511	A	G	-1	-1	6	4
chr1	69761	A	T	1	7	12	9
chr1	69871	G	A	-1	4	2	2" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_depths.1094PC0005, gt_depths.1094PC0009, \
	                    gt_depths.1094PC0012, gt_depths.1094PC0013 \
	             from variants" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 6.  Test accessing individual genotype ref_depths
####################################################################
#GT:AD:DP:GQ:PL	./.	./.	./.	./.
#GT:AD:DP:GQ:PL	0/0:7,0:7:15.04:0,15,177	0/0:2,0:2:3.01:0,3,39	0/0:6,0:6:12.02:0,12,143	0/0:4,0:4:9.03:0,9,119
#GT:AD:DP:GQ:PL	0/0:8,0:8:18.06:0,18,283	0/0:3,0:3:6.01:0,6,60	0/1:5,1:6:17.12:17,0,215	0/0:5,0:5:9.03:0,9,142
#GT:AD:DP:GQ:PL	0/1:7,1:8:8.58:9,0,188	0/1:1,2:3:29.23:33,0,29	0/0:6,0:6:18.04:0,18,206	0/0:5,0:5:12.03:0,12,139
#GT:AD:DP:GQ:PL	./.	./.	./.	./.
#GT:AD:DP:GQ:PL	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0
#GT:AD:DP:GQ:PL	0/0:2,0:2:6.01:0,6,62	0/0:79,0:79:99:0,201,1853	0/0:87,0:87:99:0,243,2236	0/0:107,0:107:99:0,294,2744
#GT:AD:DP:GQ:PL	./.	./.	0/1:2,4:6:15.70:16,0,40	0/1:2,2:4:21.59:22,0,40
#GT:AD:DP:GQ:PL	0/0:1,0:1:3.01:0,3,33	0/1:6,1:7:12.39:12,0,177	0/0:12,0:12:36.11:0,36,442	0/0:9,0:9:27.08:0,27,324
#GT:AD:DP:GQ:PL	./.	0/0:4,0:4:12.02:0,12,129	0/0:2,0:2:6.02:0,6,67	0/0:2,0:2:6.02:0,6,73
echo "    genotypes.t06...\c"
echo "chrom	end	ref	alt	gt_ref_depths.1094PC0005	gt_ref_depths.1094PC0009	gt_ref_depths.1094PC0012	gt_ref_depths.1094PC0013
chr1	30548	T	G	-1	-1	-1	-1
chr1	30860	G	C	7	2	6	4
chr1	30869	CCT	C	8	3	5	5
chr1	30895	T	C	7	1	6	5
chr1	30923	G	T	-1	-1	-1	-1
chr1	69270	A	G	-1	-1	0	0
chr1	69428	T	G	2	79	87	107
chr1	69511	A	G	-1	-1	2	2
chr1	69761	A	T	1	6	12	9
chr1	69871	G	A	-1	4	2	2" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_ref_depths.1094PC0005, gt_ref_depths.1094PC0009, \
	                    gt_ref_depths.1094PC0012, gt_ref_depths.1094PC0013 \
	             from variants" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

####################################################################
# 7.  Test accessing individual genotype alt_depths
####################################################################
#GT:AD:DP:GQ:PL	./.	./.	./.	./.
#GT:AD:DP:GQ:PL	0/0:7,0:7:15.04:0,15,177	0/0:2,0:2:3.01:0,3,39	0/0:6,0:6:12.02:0,12,143	0/0:4,0:4:9.03:0,9,119
#GT:AD:DP:GQ:PL	0/0:8,0:8:18.06:0,18,283	0/0:3,0:3:6.01:0,6,60	0/1:5,1:6:17.12:17,0,215	0/0:5,0:5:9.03:0,9,142
#GT:AD:DP:GQ:PL	0/1:7,1:8:8.58:9,0,188	0/1:1,2:3:29.23:33,0,29	0/0:6,0:6:18.04:0,18,206	0/0:5,0:5:12.03:0,12,139
#GT:AD:DP:GQ:PL	./.	./.	./.	./.
#GT:AD:DP:GQ:PL	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0
#GT:AD:DP:GQ:PL	0/0:2,0:2:6.01:0,6,62	0/0:79,0:79:99:0,201,1853	0/0:87,0:87:99:0,243,2236	0/0:107,0:107:99:0,294,2744
#GT:AD:DP:GQ:PL	./.	./.	0/1:2,4:6:15.70:16,0,40	0/1:2,2:4:21.59:22,0,40
#GT:AD:DP:GQ:PL	0/0:1,0:1:3.01:0,3,33	0/1:6,1:7:12.39:12,0,177	0/0:12,0:12:36.11:0,36,442	0/0:9,0:9:27.08:0,27,324
#GT:AD:DP:GQ:PL	./.	0/0:4,0:4:12.02:0,12,129	0/0:2,0:2:6.02:0,6,67	0/0:2,0:2:6.02:0,6,73
echo "    genotypes.t07...\c"
echo "chrom	end	ref	alt	gt_alt_depths.1094PC0005	gt_alt_depths.1094PC0009	gt_alt_depths.1094PC0012	gt_alt_depths.1094PC0013
chr1	30548	T	G	-1	-1	-1	-1
chr1	30860	G	C	0	0	0	0
chr1	30869	CCT	C	0	0	1	0
chr1	30895	T	C	1	2	0	0
chr1	30923	G	T	-1	-1	-1	-1
chr1	69270	A	G	-1	-1	3	6
chr1	69428	T	G	0	0	0	0
chr1	69511	A	G	-1	-1	4	2
chr1	69761	A	T	0	1	0	0
chr1	69871	G	A	-1	0	0	0" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_alt_depths.1094PC0005, gt_alt_depths.1094PC0009, \
	                    gt_alt_depths.1094PC0012, gt_alt_depths.1094PC0013 \
	             from variants" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp


####################################################################
# 8.  Test accessing individual genotype qualities
####################################################################
#GT:AD:DP:GQ:PL	./.	./.	./.	./.
#GT:AD:DP:GQ:PL	0/0:7,0:7:15.04:0,15,177	0/0:2,0:2:3.01:0,3,39	0/0:6,0:6:12.02:0,12,143	0/0:4,0:4:9.03:0,9,119
#GT:AD:DP:GQ:PL	0/0:8,0:8:18.06:0,18,283	0/0:3,0:3:6.01:0,6,60	0/1:5,1:6:17.12:17,0,215	0/0:5,0:5:9.03:0,9,142
#GT:AD:DP:GQ:PL	0/1:7,1:8:8.58:9,0,188	0/1:1,2:3:29.23:33,0,29	0/0:6,0:6:18.04:0,18,206	0/0:5,0:5:12.03:0,12,139
#GT:AD:DP:GQ:PL	./.	./.	./.	./.
#GT:AD:DP:GQ:PL	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0
#GT:AD:DP:GQ:PL	0/0:2,0:2:6.01:0,6,62	0/0:79,0:79:99:0,201,1853	0/0:87,0:87:99:0,243,2236	0/0:107,0:107:99:0,294,2744
#GT:AD:DP:GQ:PL	./.	./.	0/1:2,4:6:15.70:16,0,40	0/1:2,2:4:21.59:22,0,40
#GT:AD:DP:GQ:PL	0/0:1,0:1:3.01:0,3,33	0/1:6,1:7:12.39:12,0,177	0/0:12,0:12:36.11:0,36,442	0/0:9,0:9:27.08:0,27,324
#GT:AD:DP:GQ:PL	./.	0/0:4,0:4:12.02:0,12,129	0/0:2,0:2:6.02:0,6,67	0/0:2,0:2:6.02:0,6,73
echo "    genotypes.t08...\c"
echo "chrom	end	ref	alt	gt_quals.1094PC0005	gt_quals.1094PC0009	gt_quals.1094PC0012	gt_quals.1094PC0013
chr1	30548	T	G	-1.0	-1.0	-1.0	-1.0
chr1	30860	G	C	15.0399999619	3.00999999046	12.0200004578	9.02999973297
chr1	30869	CCT	C	18.0599994659	6.01000022888	17.1200008392	9.02999973297
chr1	30895	T	C	8.57999992371	29.2299995422	18.0400009155	12.029999733
chr1	30923	G	T	-1.0	-1.0	-1.0	-1.0
chr1	69270	A	G	-1.0	-1.0	9.02999973297	18.0499992371
chr1	69428	T	G	6.01000022888	99.0	99.0	99.0
chr1	69511	A	G	-1.0	-1.0	15.6999998093	21.5900001526
chr1	69761	A	T	3.00999999046	12.3900003433	36.1100006104	27.0799999237
chr1	69871	G	A	-1.0	12.0200004578	6.01999998093	6.01999998093" > exp
gemini query -q "select chrom, end, ref, alt, \
	                    gt_quals.1094PC0005, gt_quals.1094PC0009, \
	                    gt_quals.1094PC0012, gt_quals.1094PC0013 \
	             from variants" \
			 --header \
			 test.snpeff.vcf.db \
       > obs
check obs exp
rm obs exp

echo "    genotypes.t09...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" extended_ped.db --gt-filter "(gt_types).(phenotype == 1).(==HOM_REF).(all) and ((gt_types).(phenotype==2).(==HET).(any) or (gt_types).(phenotype==2).(==HOM_ALT).(any))" > obs

echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	135210790	135210791	T	C	T/T	C/C	C/C	T/T
chr10	135369531	135369532	T	C	T/T	T/C	T/C	T/T" > exp
check obs exp
rm obs exp

echo "    genotypes.t10...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(*) from variants" extended_ped.db --gt-filter "(gt_types).(phenotype == 1).(==HOM_ALT).(all) and ((gt_types).(phenotype==2).(==HET).(any) or (gt_types).(phenotype==2).(==HOM_ALT).(any))" > obs

echo "chrom	start	end	ref	alt	gts.M10475	gts.M10478	gts.M10500	gts.M128215
chr10	1142207	1142208	T	C	C/C	C/C	C/C	C/C" > exp
check obs exp


echo "    genotypes.t11...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==2) from variants limit 1" extended_ped.db > obs
echo "chrom	start	end	ref	alt	gts.M10478	gts.M10500
chr10	1142207	1142208	T	C	C/C	C/C" > exp
check obs exp


echo "    genotypes.t12...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==2), (gt_ref_depths).(phenotype==2) from variants"     --gt-filter "(gt_ref_depths).(phenotype==2).(>=20).(all)" extended_ped.db > obs

echo "chrom	start	end	ref	alt	gts.M10478	gts.M10500	gt_ref_depths.M10478	gt_ref_depths.M10500
chr10	48003991	48003992	C	T	C/T	C/T	20	24
chr16	72057434	72057435	C	T	C/C	C/C	56	67" > exp
check obs exp



echo "    genotypes.t13...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gts).(phenotype==2), (gt_ref_depths).(phenotype==2) from variants" --gt-filter "(gt_ref_depths).(phenotype==2).(>=20).(any)" extended_ped.db > obs


echo "chrom	start	end	ref	alt	gts.M10478	gts.M10500	gt_ref_depths.M10478	gt_ref_depths.M10500
chr10	48003991	48003992	C	T	C/T	C/T	20	24
chr10	126678091	126678092	G	A	G/G	G/G	11	52
chr10	135369531	135369532	T	C	T/C	T/C	27	15
chr16	72057434	72057435	C	T	C/C	C/C	56	67" > exp
check obs exp


echo "    genotypes.t14...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gt_types).(phenotype==2) from variants" extended_ped.db > obs
echo "chrom	start	end	ref	alt	gt_types.M10478	gt_types.M10500
chr10	1142207	1142208	T	C	3	3
chr10	48003991	48003992	C	T	1	1
chr10	52004314	52004315	T	C	2	3
chr10	52497528	52497529	G	C	3	3
chr10	126678091	126678092	G	A	0	0
chr10	135210790	135210791	T	C	3	3
chr10	135336655	135336656	G	A	3	2
chr10	135369531	135369532	T	C	1	1
chr16	72057434	72057435	C	T	0	0" >exp
check obs exp


echo "    genotypes.t15...\c"
gemini query --header -q "select chrom, start, end, ref, alt, (gt_types).(phenotype==2) from variants" extended_ped.db --gt-filter "(gt_types).(phenotype==2).(==HOM_REF).(all)" > obs

echo "chrom	start	end	ref	alt	gt_types.M10478	gt_types.M10500
chr10	126678091	126678092	G	A	0	0
chr16	72057434	72057435	C	T	0	0" > exp
check obs exp


n=15
for a in gt_ref_depths gt_alt_depths gt_depths gt_types; do
	for b in HOM_REF HET HOM_ALT UNKNOWN; do
		  n=$(($n + 1))
	      echo "    genotypes.t$n...\c"
		  rm -f obs
          gemini query --header -q "select chrom, start, end, ref, alt, ($a).(phenotype==2) from variants" extended_ped.db --gt-filter "(gt_types).(phenotype==2).(==$b).(all)" > obs
		  v=$(cat obs | wc -l)
		  if [ $v -gt 1 ]; then
			  echo "ok"
		  elif [ $b = UNKNOWN ]; then
			  echo "ok"
		  else
			  echo "   FAIL", $a, $b
          fi
	done
done



