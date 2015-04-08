set -e
wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI.snps_indels.vcf.tar.gz

# extract the VCF for each individual chromosome
tar -zxvf ESP6500SI.snps_indels.vcf.tar.gz 

# combine each chrom file into a single, genome wide file
(grep ^# ESP6500SI.chr1.snps_indels.vcf; cat ESP6500SI.*.vcf | \
 grep -v ^\# | sort -k1,1 -k2,2n) \
 > ESP6500SI.all.snps_indels.vcf

# compress and tabix
bgzip ESP6500SI.all.snps_indels.vcf
tabix -p vcf ESP6500SI.all.snps_indels.vcf.gz

REF=/data/hs37d5.fa
rm -f ESP6500SI.all.snps_indels.tidy.vcf.gz
# gemini version 0.12
zless ESP6500SI.all.snps_indels.vcf.gz | python sanitize-esp.py | bgzip -c > /tmp/t.gz
tabix /tmp/t.gz
vt decompose -s /tmp/t.gz \
	| vt normalize -r $REF - | perl -pe 's/([EA_|T|AA_])AC,Number=R,Type=Integer/\1AC,Number=R,Type=String/' \
	| bgzip -c > ESP6500SI.all.snps_indels.tidy.v2.vcf.gz

tabix ESP6500SI.all.snps_indels.tidy.v2.vcf.gz
rm /tmp/t.gz
