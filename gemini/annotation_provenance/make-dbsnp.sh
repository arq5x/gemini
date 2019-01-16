#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz \
#	-O dbsnp.137.vcf.gz

#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi \
#	-O dbsnp.137.vcf.gz.tbi


#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz \
#	-O dbsnp.138.vcf.gz

#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi \
#	-O dbsnp.138.vcf.gz.tbi


#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz -O dbsnp.b141.vcf.gz
#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/All.vcf.gz -O dbsnp.b141.20140813.vcf.gz
#tabix dbsnp.b141.20140813.vcf.gz

#vt decompose -s dbsnp.b141.20140813.vcf.gz | vt normalize -r /data/human/b37/human_g1k_v37_decoy.fasta - \
#	| bgzip -c > dbsnp.b141.20140813.hg19.tidy.vcf.gz
#tabix -p vcf dbsnp.b141.20140813.hg19.tidy.vcf.gz
set -euo pipefail


# 0.19.1 downloaded on 7/25/2016
#wget -O dbsnp.b147.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz
# 0.20.2 downloaded on 1/14/2018
v=151
f=dbsnp.b$v.vcf.gz
#wget -O $f ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
#wget -O $f.tbi ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi
date=20180423

fields=$(bcftools view -h $f | grep INFO | grep -Po "ID=[^,]+" | awk '{ print "INFO/"substr($0, 4, length($0)) }'| perl -pe 's/\n/,/g' | perl -pe 's/,$//')

bcftools annotate -x $fields $f --threads 3 \
    | vt decompose -s - \
	| vt normalize -r /data/human/g1k_v37_decoy.fa - \
    | bcftools annotate -x INFO/OLD_MULTIALLELIC -O z -o dbsnp.$v.$date.tidy.vcf.gz \
	&& tabix dbsnp.$v.$date.tidy.vcf.gz 

