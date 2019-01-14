# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20130118.vcf.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20130118.vcf.gz.tbi

# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20131230.vcf.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20131230.vcf.gz.tbi

# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20140303.vcf.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf/clinvar_20140303.vcf.gz.tbi

#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20140807.vcf.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20140807.vcf.gz.tbi

#wget  ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20150305.vcf.gz \
#tabix clinvar_20150305.vcf.gz

set -euo pipefail

DATE=20190102
fasta=/data/human/g1k_v37_decoy.fa

rm -rf cv.vcf.gz
wget -O cv.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2019/clinvar_${DATE}.vcf.gz

zcat cv.vcf.gz \
		   | python clinvar.py \
           | grep -v ^NW \
		   | vt normalize -n -r $fasta - \
		   | bgzip -c > clinvar_$DATE.tidy.vcf.gz
tabix clinvar_$DATE.tidy.vcf.gz
