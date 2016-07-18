check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
		exit 1
    fi
}
export -f check

echo "    x_linked_rec.t1...\c"
echo "gene	chrom	start	end	ref	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count
ASAH2C	chrX	48004991	48004992	C	T	2	1	1_dad(1_dad;unaffected;male),1_mom(1_mom;unaffected;female),1_kid(1_kid;affected;female)	C/C,C/T,T/T	1_kid	1" > exp
gemini x_linked_recessive  \
    --columns "gene, chrom, start, end, ref, alt" \
    test.x_linked.db > obs
check obs exp

echo "    x_linked_rec.t2...\c"
rm exp && touch exp
gemini x_linked_recessive  \
    --columns "gene, chrom, start, end, ref, alt" \
	-d 1000 \
    test.x_linked.db > obs
check obs exp


echo "    x_linked_rec.t3...\c"
rm exp && touch exp
gemini x_linked_recessive  \
    --columns "gene, chrom, start, end, ref, alt" \
	--min-gq 1000 \
    test.x_linked.db > obs
check obs exp

echo "    x_linked_dom.t1...\c"
rm exp; touch exp
#echo "gene	chrom	start	end	ref	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count
#	chrX	135336655	135336656	G	A	3	1	1_dad(1_dad;unaffected;male),1_mom(1_mom;unaffected;female),1_kid(1_kid;affected;female)	G/G,G/G,G/A	1_kid	1" > exp
# de novo is not x-linked dominant
gemini x_linked_dominant  \
    --columns "gene, chrom, start, end, ref, alt" \
    test.x_linked.db > obs
check obs exp


echo "    x_linked_dom.t2...\c"
rm exp && touch exp
gemini x_linked_dominant  \
    --columns "gene, chrom, start, end, ref, alt" \
	-d 1000 \
    test.x_linked.db > obs
check obs exp


echo "gene	chrom	start	end	ref	alt	variant_id	family_id	family_members	family_genotypes	samples	family_count
	chrX	135336655	135336656	G	A	3	1	1_dad(1_dad;unaffected;male),1_mom(1_mom;unaffected;female),1_kid(1_kid;affected;female)	G/G,G/G,G/A	1_kid	1" > exp
echo "    x_linked_de_novo.t1\c"
gemini x_linked_de_novo  \
    --columns "gene, chrom, start, end, ref, alt" \
    test.x_linked.db > obs
check obs exp

rm obs exp

