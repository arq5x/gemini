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
chr1	120612001	120612005	COSM674994,COSM132830
chr5	112175900	112175904	COSM33672,COSM27546,COSM1432432,COSM25817,COSM1432432,COSM25817,COSM19698,COSM177897
chr9	139277993	139278022	COSM247490" > exp

gemini query -q "select chrom, start, end, cosmic_ids from variants" test.cosmic.db > obs
check obs exp
rm obs exp
