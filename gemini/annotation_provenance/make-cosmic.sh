# this uses CloudBioLinux's pre-cleaned and merged COSMIC VCF:
# https://github.com/chapmanb/cloudbiolinux/blob/master/utils/prepare_cosmic.py
#wget https://s3.amazonaws.com/biodata/variants/cosmic-v67_20131024-hg19.vcf.gz

#mv cosmic-v67_20131024-hg19.vcf.gz hg19.cosmic.v67.20131024.gz
#tabix -p hg19.cosmic.v67.20131024.gz

# Version 68
wget https://s3.amazonaws.com/biodata/variants/cosmic-v68-GRCh37.vcf.gz
#wget https://s3.amazonaws.com/biodata/variants/cosmic-v68-GRCh37.vcf.gz.tbi

# gemini version 0.12
vt decompose -s cosmic-v68-GRCh37.vcf.gz \
	| vt normalize -r /data/human/b37/human_g1k_v37_decoy.fasta - \
	| bgzip -c > cosmic-v68-GRCh37.tidy.vcf.gz
tabix cosmic-v68-GRCh37.tidy.vcf.gz
