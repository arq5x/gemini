#curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz \
#	> dbsnp.137.vcf.gz

#curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi \
#	> dbsnp.137.vcf.gz.tbi


# curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz \
#	> dbsnp.138.vcf.gz

# curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi \
#	> dbsnp.138.vcf.gz.tbi


curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz > dbsnp.b141.vcf.gz
gunzip dbsnp.b141.vcf.gz | bgzip -c > dbsnp.hg19.b141.vcf.gz
tabix -p vcf dbsnp.hg19.b141.vcf.gz