check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

echo " nogls.t1"
python drop-gl-columns.py test.query.db test.nogls.db
gemini bcolz_index test.nogls.db
sed 's/test.query.db/test.nogls.db/g' test-query.sh | bash | sed 's/query/nogls-sub/' 
echo "NOTE: we expect t14 to fail due to LIMIT"
sed 's/test.query.db/test.nogls.db/g' test-query.sh | sed 's/gemini query/gemini query --use-bcolz/' | bash | sed 's/query/nogls-sub-bcolz/'
wait

echo " nogls.t2"
# mendel error uses gl columns to calculate violation prob:
python drop-gl-columns.py test.mendel.db test.mendelnogls.db
rm -f obs exp
echo "chrom	start	end	variant_id	family_id	family_members	family_genotypes	samples	family_count	violation
chr1	10670	10671	1	CEPH1463	NA12889(NA12889;unknown;male),NA12890(NA12890;unknown;female),NA12877(NA12877;unknown;male)	G/G,G/G,G/C	NA12877	1	plausible de novo
chr1	28493	28494	2	CEPH1463	NA12889(NA12889;unknown;male),NA12890(NA12890;unknown;female),NA12877(NA12877;unknown;male)	T/C,T/T,C/C	NA12877	1	loss of heterozygosity
chr1	28627	28628	3	CEPH1463	NA12889(NA12889;unknown;male),NA12890(NA12890;unknown;female),NA12877(NA12877;unknown;male)	C/C,C/C,C/T	NA12877	1	plausible de novo" > exp
gemini mendel_errors --columns "chrom,start,end" test.mendelnogls.db | head -4 > obs
check obs exp
rm obs exp

echo "nogls.comphet.t3"

python drop-gl-columns.py test.comp_het_default.db test.comp_het.nogls.db
echo "chrom	start	end	ref	alt	gene	impact	variant_id	family_id	family_members	family_genotypes	samples	family_count	comp_het_id
chr1	17362	17366	TTCT	T	WASH7P	splice_acceptor	1	0	child_1(child_1;unknown),child_2(child_2;unknown),dad_2(dad_2;unknown),mom_2(mom_2;unknown),dad_1(dad_1;unknown),mom_1(mom_1;unknown),child_3(child_3;unknown),dad_3(dad_3;unknown),mom_3(mom_3;unknown),child_4(child_4;unknown),dad_4(dad_4;unknown),mom_4(mom_4;unknown)	TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|TTCT,TTCT|T,TTCT|T,TTCT|TTCT		1	1_1_2
chr1	17729	17730	C	A	WASH7P	splice_acceptor	2	0	child_1(child_1;unknown),child_2(child_2;unknown),dad_2(dad_2;unknown),mom_2(mom_2;unknown),dad_1(dad_1;unknown),mom_1(mom_1;unknown),child_3(child_3;unknown),dad_3(dad_3;unknown),mom_3(mom_3;unknown),child_4(child_4;unknown),dad_4(dad_4;unknown),mom_4(mom_4;unknown)	C|C,C|C,C|C,C|C,C|A,C|C,C|A,C|C,C|C,A|C,C|A,C|A		1	1_1_2" > exp

gemini comp_hets \
    --column "chrom,start,end,ref,alt,gene,impact" \
    --allow-unaffected \
    --lenient \
    test.comp_het.nogls.db > obs
check obs exp
rm obs exp


echo "nogls.auto_dom.t4"
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	3	3	3_dad(3_dad;affected),3_mom(3_mom;unknown),3_kid(3_kid;affected)	C/T,C/C,C/T	3_dad,3_kid	2
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	3	2	2_dad(2_dad;unaffected),2_mom(2_mom;affected),2_kid(2_kid;affected)	C/C,C/T,C/T	2_mom,2_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	4	3	3_dad(3_dad;affected),3_mom(3_mom;unknown),3_kid(3_kid;affected)	C/T,C/C,C/T	3_dad,3_kid	2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	4	2	2_dad(2_dad;unaffected),2_mom(2_mom;affected),2_kid(2_kid;affected)	C/C,C/T,C/T	2_mom,2_kid	2
SPRN	chr10	135336655	135336656	G	A	intron	LOW	5	3	3_dad(3_dad;affected),3_mom(3_mom;unknown),3_kid(3_kid;affected)	G/A,G/G,G/A	3_dad,3_kid	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	3	3_dad(3_dad;affected),3_mom(3_mom;unknown),3_kid(3_kid;affected)	T/C,T/T,T/C	3_dad,3_kid	2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	2_dad(2_dad;unaffected),2_mom(2_mom;affected),2_kid(2_kid;affected)	T/T,T/C,T/C	2_mom,2_kid	2
WDR37	chr10	1142208	1142209	T	C	stop_loss	HIGH	2	3	3_dad(3_dad;affected),3_mom(3_mom;unknown),3_kid(3_kid;affected)	T/C,T/C,T/C	3_dad,3_kid	2" > exp

python drop-gl-columns.py test.auto_dom.db test.aa.db
gemini autosomal_dominant  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
    --lenient \
    test.aa.db > obs
check obs exp
rm obs exp test.aa.db

