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
chr1	10670	10671	1	CEPH1463	NA12889(dad;unknown),NA12890(mom;unknown),NA12877(child;unknown)	G/G,G/G,G/C	NA12877	1	plausible de novo	0.962
chr1	28493	28494	2	CEPH1463	NA12889(dad;unknown),NA12890(mom;unknown),NA12877(child;unknown)	T/C,T/T,C/C	NA12877	1	loss of heterozygosity	0.660
chr1	28627	28628	3	CEPH1463	NA12889(dad;unknown),NA12890(mom;unknown),NA12877(child;unknown)	C/C,C/C,C/T	NA12877	1	plausible de novo	0.989" > exp

check obs exp


###############
# 2. Test depth
###############

gemini mendel_errors --columns "chrom,start,end" test.mendel.db -d 1000 > obs
echo "chrom	start	end	variant_id	family_id	family_members	family_genotypes	samples	family_count	violation	violation_prob" > exp
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

