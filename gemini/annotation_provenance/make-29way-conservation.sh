#/usr/bin/sh

# grab the 29-way conservation file from Broad
wget http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz
# conservation file is in hg18. liftover to hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
liftover 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz \
         hg18ToHg19.over.chain.gz \
         29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.bed \
         unmapped
# merge overlapping intervals
sort -k1,1 -k2,2n 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.bed | \
   bedtools merge -i  - \
   > 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed
rm 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.bed
# gbzip for Tabix
bgzip 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed
# tabix
tabix -p bed 29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.hg19.merged.bed.gz
