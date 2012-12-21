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