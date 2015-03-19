#############################################
# 1. Test ExAC allele frequencies
#############################################
check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check


echo "    exac.t01...\c"
echo "chrom	start	end	ref	alt	in_exac	aaf_exac_all	aaf_adj_exac_all	aaf_adj_exac_afr	aaf_adj_exac_amr	aaf_adj_exac_eas	aaf_adj_exac_fin	aaf_adj_exac_nfe	aaf_adj_exac_oth	aaf_adj_exac_sas
chr1	985954	985955	G	C	0	None	None	None	None	None	None	None	None	None
chr1	1199488	1199489	G	A	0	None	None	None	None	None	None	None	None	None
chr1	1959698	1959699	G	A	1	0.017	0.0175763511709	0.00489987217725	0.00975699558174	0	0.0412407472682	0.0237794612795	0.0224438902743	0.00647021140295
chr1	161276552	161276553	G	T	0	None	None	None	None	None	None	None	None	None
chr1	247587092	247587093	C	T	0	None	None	None	None	None	None	None	None	None
chr9	112185055	112185056	C	G	1	8.236e-06	8.23709658819e-06	0	0	0	0	1.49853143919e-05	0	0
chr12	121432116	121432118	GC	G	1	0.243	0.278784178031	0.215484048405	0.486200895077	0.00100267379679	0.317728531856	0.330323021843	0.255813953488	0.109512269009
chr14	105420589	105420590	C	T	1	0.001735	0.00174129353234	0.00030637254902	0.000777873811582	0	0.00030238887209	0.00246261036699	0.00447427293065	0.00169594185342" > exp

gemini query --header -q "select chrom, start, end, ref, alt, in_exac, aaf_exac_all, aaf_adj_exac_all, \
	                        aaf_adj_exac_afr, aaf_adj_exac_amr, aaf_adj_exac_eas, aaf_adj_exac_fin, \
		                        aaf_adj_exac_nfe, aaf_adj_exac_oth, aaf_adj_exac_sas from variants" test.exac.db > obs


check obs exp
rm obs exp

#############################################
