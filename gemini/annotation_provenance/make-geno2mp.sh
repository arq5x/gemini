# 0.18.1
wget -O - http://geno2mp.gs.washington.edu/download/Geno2MP.variants.vcf \
	 | vt decompose -s - \
	 | vt normalize -r /data/human/b37/human_g1k_v37_decoy.fasta - \
	 | bgzip -c > geno2mp.variants.tidy.vcf.gz
tabix geno2mp.variants.tidy.vcf.gz
