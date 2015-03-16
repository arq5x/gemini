#curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz \
#	> dbsnp.137.vcf.gz

#curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi \
#	> dbsnp.137.vcf.gz.tbi


# curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz \
#	> dbsnp.138.vcf.gz

# curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi \
#	> dbsnp.138.vcf.gz.tbi


# curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz > dbsnp.b141.vcf.gz
#curl ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/All.vcf.gz > dbsnp.b141.20140813.vcf.gz
tabix dbsnp.b141.20140813.vcf.gz

vt decompose -s dbsnp.b141.20140813.vcf.gz | vt normalize -r /data/human/b37/human_g1k_v37_decoy.fasta - \
	| bgzip -c > dbsnp.b141.20140813.hg19.tidy.vcf.gz
tabix -p vcf dbsnp.b141.20140813.hg19.tidy.vcf.gz
