#######################################################################################
# 1.Test gemini region (--reg)
#######################################################################################
echo "    region.t01...\c"
echo "chr1	10000	10001	T	TC	DDX11L1
chr1	10055	10056	A	C	DDX11L1" > exp

gemini region --reg chr1:10000-10100 --columns "chrom, start, end, ref, alt, gene"  test.region.db > obs

check obs exp
rm obs exp

#######################################################################################
# 2.Test gemini region (--columns)
#######################################################################################
echo "    region.t02...\c"
echo "chr16	72057281	72057282	A	G	DHODH
chr16	72057434	72057435	C	T	DHODH
chr16	72059268	72059269	T	C	DHODH" > exp

gemini region --gene DHODH --columns "chrom, start, end, ref, alt, gene" test.region.db > obs
check obs exp
rm obs exp

#######################################################################################
# 3.Test gemini region (--columns and --filter)
#######################################################################################
echo "    region.t03...\c"
echo "chr16	72057281	72057282	A	G	DHODH" > exp

gemini region --gene DHODH --columns "chrom, start, end, ref, alt, gene" --filter "alt='G'" test.region.db > obs
check obs exp
rm obs exp

#######################################################################################
# 4. Test gemini region (--columns and --filter and --header)
#######################################################################################
echo "    region.t04...\c"
echo "chrom	start	end	ref	alt	gene
chr16	72057281	72057282	A	G	DHODH" > exp

gemini region --gene DHODH --columns "chrom, start, end, ref, alt, gene" --filter "alt='G'" --header test.region.db > obs
check obs exp
rm obs exp

#######################################################################################
# 5. Test gemini region (--columns and --filter and --json)
#######################################################################################
echo "    region.t05...\c"
echo "{\"end\": 72057282, \"ref\": \"A\", \"start\": 72057281, \"alt\": \"G\", \"gene\": \"DHODH\", \"chrom\": \"chr16\"}" > exp

gemini region --format json --gene DHODH --columns "chrom, start, end, ref, alt, gene" --filter "alt='G'" test.region.db > obs
check obs exp
rm obs exp
