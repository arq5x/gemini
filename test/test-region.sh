#######################################################################################
# 1.Test gemini region (--reg)
#######################################################################################
echo "    region.t01...\c"
echo "chrom,start,end,variant_id,anno_id,ref,alt,qual,filter,type,sub_type,\
call_rate,in_dbsnp,rs_ids,in_omim,clinvar_sig,clinvar_disease_name,\
clinvar_dbsource,clinvar_dbsource_id,clinvar_origin,clinvar_dsdb,\
clinvar_dsdbid,clinvar_disease_acc,clinvar_in_locus_spec_db,\
clinvar_on_diag_assay,pfam_domain,cyto_band,rmsk,in_cpg_island,in_segdup,is_conserved,\
gerp_bp_score,gerp_element_pval,num_hom_ref,num_het,num_hom_alt,num_unknown,aaf,hwe,\
inbreeding_coeff,pi,recomb_rate,gene,transcript,is_exonic,is_coding,is_lof,exon,\
codon_change,aa_change,aa_length,biotype,impact,impact_severity,polyphen_pred,\
polyphen_score,sift_pred,sift_score,anc_allele,rms_bq,cigar,depth,strand_bias,\
rms_map_qual,in_hom_run,num_mapq_zero,num_alleles,num_reads_w_dels,\
haplotype_score,qual_depth,allele_count,allele_bal,in_hm2,in_hm3,is_somatic,\
in_esp,aaf_esp_ea,aaf_esp_aa,aaf_esp_all,exome_chip,in_1kg,aaf_1kg_amr,\
aaf_1kg_asn,aaf_1kg_afr,aaf_1kg_eur,aaf_1kg_all,grc,gms_illumina,gms_solid,\
gms_iontorrent,in_cse,encode_tfbs,encode_dnaseI_cell_count,encode_dnaseI_cell_list,\
encode_consensus_gm12878,encode_consensus_h1hesc,encode_consensus_helas3,\
encode_consensus_hepg2,encode_consensus_huvec,encode_consensus_k562
chr1,10000,10001,1,1,T,TC,175.91,.,indel,ins,1.0,0,.,.,.,.,.,.,.,.,.,.,.,.,.,\
chr1p36.33,Simple_repeat_Simple_repeat_(CCCTAA)n;trf;Satellite_telo_TAR1;trf;trf;\
trf;trf;trf,0,1,0,.,.,0,4,0,0,0.5,0.045500265919,-1,0.571428571429,2.981822,DDX11L1,\
ENST00000456328,0,0,0,.,.,.,.,.,upstream,LOW,.,.,.,.,.,.,.,76,.,35.31,0,0,8,.,\
218.6157,2.31,4,.,.,.,.,0,.,.,.,0,0,.,.,.,.,.,.,.,.,.,0,.,.,.,CTCF,CTCF,unknown,\
unknown,unknown,CTCF
chr1,10055,10056,2,1,A,C,47.27,.,snp,tv,1.0,0,.,.,.,.,.,.,.,.,.,.,.,.,.,chr1p36.33,\
Simple_repeat_Simple_repeat_(CCCTAA)n;trf;Satellite_telo_TAR1;trf;trf;trf;trf;trf,\
0,1,0,.,.,3,1,0,0,0.125,0.775096962148,-0.142857142857,0.25,2.981822,DDX11L1,\
ENST00000456328,0,0,0,.,.,.,.,.,upstream,LOW,.,.,.,.,.,.,.,199,.,30.59,0,0,8,0.0,\
10.0031,0.65,1,.,.,.,.,0,.,.,.,0,0,.,.,.,.,.,.,89.8,29.5,92.8,0,.,.,.,CTCF,CTCF,\
unknown,unknown,unknown,CTCF" > exp

gemini region --reg chr1:10000-10100 --header --sep "," test.region.db > obs

check obs exp
rm obs exp

#######################################################################################
# 2.Test gemini region (--gene)
#######################################################################################
echo "    region.t02...\c"
echo "chrom,start,end,variant_id,anno_id,ref,alt,qual,filter,type,sub_type,\
call_rate,in_dbsnp,rs_ids,in_omim,clinvar_sig,clinvar_disease_name,clinvar_dbsource,\
clinvar_dbsource_id,clinvar_origin,clinvar_dsdb,clinvar_dsdbid,clinvar_disease_acc,\
clinvar_in_locus_spec_db,clinvar_on_diag_assay,pfam_domain,cyto_band,rmsk,in_cpg_island,\
in_segdup,is_conserved,gerp_bp_score,gerp_element_pval,num_hom_ref,num_het,num_hom_alt,\
num_unknown,aaf,hwe,inbreeding_coeff,pi,recomb_rate,gene,transcript,is_exonic,is_coding,\
is_lof,exon,codon_change,aa_change,aa_length,biotype,impact,impact_severity,polyphen_pred,\
polyphen_score,sift_pred,sift_score,anc_allele,rms_bq,cigar,depth,strand_bias,rms_map_qual,\
in_hom_run,num_mapq_zero,num_alleles,num_reads_w_dels,haplotype_score,qual_depth,allele_count,\
allele_bal,in_hm2,in_hm3,is_somatic,in_esp,aaf_esp_ea,aaf_esp_aa,aaf_esp_all,exome_chip,\
in_1kg,aaf_1kg_amr,aaf_1kg_asn,aaf_1kg_afr,aaf_1kg_eur,aaf_1kg_all,grc,gms_illumina,\
gms_solid,gms_iontorrent,in_cse,encode_tfbs,encode_dnaseI_cell_count,encode_dnaseI_cell_list,\
encode_consensus_gm12878,encode_consensus_h1hesc,encode_consensus_helas3,\
encode_consensus_hepg2,encode_consensus_huvec,encode_consensus_k562
chr16,72057281,72057282,7,1,A,G,1098.89,.,snp,ts,1.0,1,rs2288002,.,.,.,.,.,.,.,.,.,.,.,\
.,chr16q22.2,.,0,0,0,.,.,0,3,1,0,0.625,0.230139340577,-0.6,0.535714285714,0.174767,DHODH,\
ENST00000219240,0,0,0,.,.,.,.,.,intron,LOW,.,.,.,.,.,.,.,80,.,36.24,3,0,8,0.0,1.9993,\
13.74,5,.,.,.,.,0,.,.,.,0,1,0.43,0.38,0.61,0.58,0.51,.,.,.,.,0,.,.,.,T,unknown,T,T,T,T
chr16,72057434,72057435,8,1,C,T,572.98,.,snp,ts,1.0,1,rs201947120,1,pathogenic,\
Miller syndrome,OMIM Allelic Variant,126064.0001,.,NCBI:OMIM:Orphanet:SNOMED CT,\
C0265257:263750:246:66038001,RCV000018291.1,1,0,DHO_dh,chr16q22.2,.,0,0,1,.,5.8541e-64,3,1,0,\
0,0.125,0.775096962148,-0.142857142857,0.25,0.140436,DHODH,ENST00000219240,1,1,0,8/9,Cgg/Tgg,\
R/W,.,.,non_syn_coding,MED,probably_damaging,0.984,deleterious,0.0,.,.,.,260,.,36.53,\
0,0,8,0.0,4.5319,8.07,1,.,.,.,.,1,0.000239,0,0.000161,1,0,.,.,.,.,.,.,.,.,.,0,.,3,\
8988t;CACO2;T47d,T,T,T,T,T,T
chr16,72059268,72059269,9,1,T,C,39.18,.,snp,ts,0.25,1,rs2287998,.,.,.,.,.,.,.,.,.,.,.,\
.,chr16q22.2,.,0,0,0,.,.,0,0,1,3,1.0,1,.,0,0.102,DHODH,ENST00000219240,0,0,0,.,.,.,.,.,\
downstream,LOW,.,.,.,.,.,.,.,2,.,37.0,0,0,2,0.0,0.0,19.59,2,.,.,.,.,0,.,.,.,0,1,0.43,\
0.38,0.61,0.58,0.51,.,.,.,.,0,JunD_1,.,.,T,T,T,unknown,T,T" > exp

gemini region --gene DHODH --header --sep "," test.region.db > obs
check obs exp
rm obs exp