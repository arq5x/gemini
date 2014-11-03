#############################################
# 1. Test ExAC allele frequencies
#############################################
echo "    exac.t01...\c"
echo "chrom	start	end	in_exac	aaf_exac_all	aaf_adj_exac_all	aaf_adj_exac_afr	aaf_adj_exac_amr	aaf_adj_exac_eas	aaf_adj_exac_fin	aaf_adj_exac_nfe	aaf_adj_exac_oth	aaf_adj_exac_sas
chr1	985954	985955	0	None	None	None	None	None	None	None	None	None
chr1	1199488	1199489	0	None	None	None	None	None	None	None	None	None
chr1	1959698	1959699	1	0.017	0.0175534349801	0.00481877226063	0.00973012667523	0	0.0410552061495	0.0237265929571	0.0232843137255	0.00655049605698
chr1	161276552	161276553	0	None	None	None	None	None	None	None	None	None
chr1	235891374	235891431	0	None	None	None	None	None	None	None	None	None
chr1	247587092	247587093	0	None	None	None	None	None	None	None	None	None
chr3	33063140	33063141	0	None	None	None	None	None	None	None	None	None
chrX	31200831	31200832	0	None	None	None	None	None	None	None	None	None
chrX	89963313	89963314	0	None	None	None	None	None	None	None	None	None" > exp

gemini query --header -q "select chrom, start, end, in_exac, aaf_exac_all, aaf_adj_exac_all, \
	                        aaf_adj_exac_afr, aaf_adj_exac_amr, aaf_adj_exac_eas, aaf_adj_exac_fin, \
		                        aaf_adj_exac_nfe, aaf_adj_exac_oth, aaf_adj_exac_sas from variants" test.clinvar.db > obs

check obs exp
rm obs exp

#############################################
