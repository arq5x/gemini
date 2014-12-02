#############################################
# 1. Test ExAC allele frequencies
#############################################
echo "    exac.t01...\c"
echo "chrom	start	end	ref	alt	in_exac	aaf_exac_all	aaf_adj_exac_all	aaf_adj_exac_afr	aaf_adj_exac_amr	aaf_adj_exac_eas	aaf_adj_exac_fin	aaf_adj_exac_nfe	aaf_adj_exac_oth	aaf_adj_exac_sas
chr1	985954	985955	G	C	0	None	None	None	None	None	None	None	None	None
chr1	1199488	1199489	G	A	0	None	None	None	None	None	None	None	None	None
chr1	1959698	1959699	G	A	1	0.017	0.0175534349801	0.00481877226063	0.00973012667523	0	0.0410552061495	0.0237265929571	0.0232843137255	0.00655049605698
chr1	161276552	161276553	G	T	0	None	None	None	None	None	None	None	None	None
chr1	247587092	247587093	C	T	0	None	None	None	None	None	None	None	None	None
chr9	112185055	112185056	C	G	1	1.626e-05	1.62665105082e-05	0.000189000189	0	0	0	0	0	0
chr12	121432116	121432118	GC	G	1	0.0002284	0.000234270414993	0.000539374325782	0.000371839365394	0.000825354902608	0	9.89472017731e-05	0	0.000300571085062
chr14	105420589	105420590	C	T	1	0.001746	0.00175203039036	0.000301083902047	0.000775862068966	0	0.000296384113811	0.0024867521241	0.00437636761488	0.00168431183831" > exp

gemini query --header -q "select chrom, start, end, ref, alt, in_exac, aaf_exac_all, aaf_adj_exac_all, \
	                        aaf_adj_exac_afr, aaf_adj_exac_amr, aaf_adj_exac_eas, aaf_adj_exac_fin, \
		                        aaf_adj_exac_nfe, aaf_adj_exac_oth, aaf_adj_exac_sas from variants" test.exac.db > obs

check obs exp
rm obs exp

#############################################