###########################################################################################
#1. Test loading an unannotated file without genotypes
###########################################################################################
gemini load -v ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.vcf \
	        --skip-gene-tables --skip-gerp-bp --skip-cadd --no-genotypes 1000G.snippet.db

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
    --skip-gene-tables --skip-gerp-bp --skip-cadd --no-genotypes \
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

###########################################################################################
#3. Test loading an extended ped file
###########################################################################################
gemini load -p test_extended_ped.ped -v test4.vep.snpeff.vcf \
--skip-gene-tables --skip-gerp-bp --skip-cadd -t snpEff extended_ped_test.db

echo "    load.t3...\c"
echo "sample_id	family_id	name	paternal_id	maternal_id	sex	phenotype	ethnicity	hair_color
1	1	M10475	0	0	1	1	None	brown
2	1	M10478	M10475	M10500	2	2	None	brown
3	1	M10500	0	0	2	2	None	purple
4	1	M128215	M10475	M10500	1	1	None	blue" > exp
gemini query --header -q "select * from samples" extended_ped_test.db > obs
check obs exp
rm obs exp

###########################################################################################
#4. Test --passonly on loading
###########################################################################################
gemini load  --skip-gene-tables --passonly -v test.passonly.vcf --skip-gerp-bp --skip-cadd -t snpEff \
passonly.db

echo "    load.t4...\c"
echo "chr1	1334051	CTAGAG	C" > exp
gemini query -q "select chrom, start, ref, alt from variants" passonly.db > obs
check obs exp
rm obs exp
