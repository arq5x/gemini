# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20130118.vcf.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20130118.vcf.gz.tbi

# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20131230.vcf.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20131230.vcf.gz.tbi

# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20140303.vcf.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20140303.vcf.gz.tbi

#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20140807.vcf.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20140807.vcf.gz.tbi

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20150305.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20150305.vcf.gz.tbi

VCF=clinvar_20150305.vcf.gz
zless $VCF | perl -pe 's/##INFO=<ID=CLN(DBN|ACC|SRCID|SIG|DADB|SRC|ORIGIN|ALLE|HGVS|DSDBID|DSDBD|DSDB),Number=\./##INFO=<ID=CLN\1,Number=A/' | bgzip -c > tmp.vcf.gz

tabix tmp.vcf.gz

vt decompose -s tmp.vcf.gz \
		   | vt normalize -r /data/human/b37/human_g1k_v37_decoy.fasta - \
		   | bgzip -c > clinvar_20150305.tidy.vcf.gz
tabix clinvar_20150305.tidy.vcf.gz
rm tmp.vcf.gz
