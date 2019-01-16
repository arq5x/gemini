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


echo "    gnomad.t01...\n"
echo "chrom	start	end	ref	alt	aaf_gnomad_all	aaf_gnomad_all	aaf_gnomad_afr	aaf_gnomad_amr	aaf_gnomad_eas	aaf_gnomad_fin	aaf_gnomad_nfe	aaf_gnomad_oth	aaf_gnomad_sas
chr1	985954	985955	G	C	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr1	1199488	1199489	G	A	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr1	1959698	1959699	G	A	0.0172546003014	0.0172546003014	0.00383117698782	0.00964331358197	0.000109206071858	0.0381636652048	0.0233894058697	0.0202702702703	0.0061303435615
chr1	161276552	161276553	G	T	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr1	247587092	247587093	C	T	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr9	112185055	112185056	C	G	7.95558025857e-06	7.95558025857e-06	0.000123061777012	0.0	0.0	0.0	0.0	0.0	0.0
chr12	121432116	121432118	GC	G	6.13544980297e-05	6.13544980297e-05	0.000145095763204	3.11526479751e-05	5.88166098106e-05	5.12400081984e-05	7.85792865001e-05	0.000177179305457	0.0
chr14	105420589	105420590	C	T	0.00171821005642	0.00171821005642	0.000258297817383	0.000898030127462	0.0	0.000190566936636	0.00195814357357	0.00248426631335	0.00163409373162" > exp.$$

gemini query --header -q "select chrom, start, end, ref, alt, aaf_gnomad_all, aaf_gnomad_all, \
	                        aaf_gnomad_afr, aaf_gnomad_amr, aaf_gnomad_eas, aaf_gnomad_fin, \
		                        aaf_gnomad_nfe, aaf_gnomad_oth, aaf_gnomad_sas from variants" test.exac.db > obs.$$


check obs.$$ exp.$$
rm obs.$$ exp.$$

echo "    gnomad.t02...\n"
gemini query --header -q "select chrom, start, end, ref, alt, aaf_gnomad_non_neuro, aaf_gnomad_controls, aaf_gnomad_non_cancer from variants" test.exac.db > obs.$$
echo "chrom	start	end	ref	alt	aaf_gnomad_non_neuro	aaf_gnomad_controls	aaf_gnomad_non_cancer
chr1	985954	985955	G	C	-1.0	-1.0	-1.0
chr1	1199488	1199489	G	A	-1.0	-1.0	-1.0
chr1	1959698	1959699	G	A	0.0164560005069	0.0170828998089	0.0171958003193
chr1	161276552	161276553	G	T	-1.0	-1.0	-1.0
chr1	247587092	247587093	C	T	-1.0	-1.0	-1.0
chr9	112185055	112185056	C	G	9.61224031926e-06	0.0	4.22193988925e-06
chr12	121432116	121432118	GC	G	5.28017008037e-05	6.08346999798e-05	5.99675004196e-05
chr14	105420589	105420590	C	T	0.0016407299554	0.00138945004437	0.00172107003164" > exp.$$

check obs.$$ exp.$$
rm obs.$$ exp.$$

#############################################

#############################################
