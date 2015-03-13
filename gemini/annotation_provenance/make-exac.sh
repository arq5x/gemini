# v0.12
wget -O - ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz \
	| vt decompose - | vt normalize -r /data/human/b37/human_g1k_v37_decoy.fasta - \
	| bgzip -c > ExAC.r0.3.sites.vep.tidy.vcf.gz
tabix ExAC.r0.3.sites.vep.tidy.vcf.gz
