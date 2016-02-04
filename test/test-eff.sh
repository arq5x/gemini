echo "test.eff-plus.t1"
echo "chrom	start	end	gene	ref	alt	aa_change	is_coding	is_exonic	is_splicing
chr6	34950530	34950531	ANKS1A	G	A	T245	0	0	1" > exp
gemini query --header -q "select chrom,start,end,gene,ref,alt,aa_change,is_coding,is_exonic,is_splicing from variants" test.eff.db > obs

check obs exp
