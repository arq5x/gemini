check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}
export -f check
###################################################################
# 1. Test basic functionality
###################################################################
echo "    mendel_error.t1..."

gemini mendel_errors --columns "chrom,start,end" test.mendel.db | head -4 > obs
echo "chrom	start	end	variant_id	family_id	family_members	family_genotypes	samples	family_count	violation	violation_prob
chr1	10670	10671	1	CEPH1463	NA12889(NA12889;unknown;male),NA12890(NA12890;unknown;female),NA12877(NA12877;unknown;male)	G/G,G/G,G/C	NA12877	1	plausible de novo	0.96228
chr1	28493	28494	2	CEPH1463	NA12889(NA12889;unknown;male),NA12890(NA12890;unknown;female),NA12877(NA12877;unknown;male)	T/C,T/T,C/C	NA12877	1	loss of heterozygosity	0.65973
chr1	28627	28628	3	CEPH1463	NA12889(NA12889;unknown;male),NA12890(NA12890;unknown;female),NA12877(NA12877;unknown;male)	C/C,C/C,C/T	NA12877	1	plausible de novo	0.98859" > exp

check obs exp
rm exp obs


###############
# 2. Test depth
###############

gemini mendel_errors --columns "chrom,start,end" test.mendel.db -d 1000 > obs
touch exp
echo "    mendel_error.t2..."
check obs exp

################
# 3-4 Test gt-ll
################

gemini mendel_errors --columns "chrom,start,end" test.mendel.db --gt-pl-max 1 | wc -l > obs
echo "11" > exp
echo "    mendel_error.t3..."
check obs exp

# with a high enough limit, we get all the variants
gemini mendel_errors --columns "chrom,start,end" test.mendel.db --gt-pl-max 1000 | wc -l > obs
echo "22" > exp
echo "    mendel_error.t4..."
check obs exp

###############
# test families
###############

gemini mendel_errors --columns "chrom,start,end" test.mendel.db --families CEPH1463 | wc -l > obs
echo "22" > exp
echo "    mendel_error.t5..."
check obs exp

