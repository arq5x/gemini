####################################################################
# 1. Test cosmic annotations
####################################################################
check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

echo "    cosmic.t01...\c"
echo "chr1	10378	10379	None
chr1	120612001	120612004	COSM132830
chr5	112175900	112175903	COSM19698
chr9	139277993	139277997	COSM247490" > exp

gemini query -q "select chrom, start, end, cosmic_ids from variants" test.cosmic.db > obs
check obs exp
rm obs exp
