#!/usr/bin/sh
# Note, despite the date in the filename (20101123), the last modified
# timestamp on the 1000G site for this download was 10/12/12 1:28:00 PM 
#curl  ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz \
#> ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.2012Oct12.vcf.gz

#curl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz.tbi > ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.2012Oct12.vcf.gz.tbi

# 1000g phase 3 data
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz > ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi > ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz.tbi

