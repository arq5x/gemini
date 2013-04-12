###########################################################################################
#1. Test loading an unannotated file without genotypes
###########################################################################################
gemini load -v ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.vcf \
	        --no-genotypes 1000G.snippet.db

echo "    load.t1...\c"
echo "chr1	10582	G	A	None
chr1	10610	C	G	None
chr1	13301	C	T	None
chr1	13326	G	C	None
chr1	13956	TC	T	None
chr1	13979	T	C	None
chr1	30922	G	T	None
chr1	46401	C	CTGT	None
chr1	47189	G	GA	None
chr1	51475	T	C	None" > exp

gemini query -q "select chrom, start, ref, alt, gene from variants limit 10" \
	1000G.snippet.db > obs
check obs exp
rm obs exp


###########################################################################################
#2. Test loading an annotated file without genotypes
###########################################################################################
gemini load -v ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.snpEff.vcf \
    --no-genotypes \
    -t snpEff 1000G.snippet.snpEff.db

echo "    load.t2...\c"
echo "chr1	10582	G	A	WASH7P
chr1	10610	C	G	WASH7P
chr1	13301	C	T	WASH7P
chr1	13326	G	C	WASH7P
chr1	13956	TC	T	DDX11L1
chr1	13979	T	C	DDX11L1
chr1	30922	G	T	FAM138A
chr1	46401	C	CTGT	None
chr1	47189	G	GA	None
chr1	51475	T	C	None" > exp

gemini query -q "select chrom, start, ref, alt, gene from variants limit 10" \
	1000G.snippet.snpEff.db > obs
check obs exp
rm obs exp