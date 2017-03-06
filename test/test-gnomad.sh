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


echo "    gnomad.t01...\c"
echo "chrom	start	end	ref	alt	aaf_gnomad_all	aaf_gnomad_all	aaf_gnomad_afr	aaf_gnomad_amr	aaf_gnomad_eas	aaf_gnomad_fin	aaf_gnomad_nfe	aaf_gnomad_oth	aaf_gnomad_sas
chr1	985954	985955	G	C	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr1	1199488	1199489	G	A	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr1	1959698	1959699	G	A	0.0174018	0.0174018	0.0036564286664	0.00955851604038	0.000116387337058	0.0382032667877	0.0233466915698	0.0213392200147	0.00599856556041
chr1	161276552	161276553	G	T	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr1	247587092	247587093	C	T	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0	-1.0
chr9	112185055	112185056	C	G	8.12315e-06	8.12315e-06	0.000130701869037	0.0	0.0	0.0	0.0	0.0	0.0
chr12	121432116	121432118	GC	G	6.27201e-05	6.27201e-05	0.000154202004626	3.19570497252e-05	6.25312656328e-05	5.06790999392e-05	8.0064051241e-05	0.000198098256735	0.0
chr14	105420589	105420590	C	T	0.00166166	0.00166166	0.000261574679571	0.000863609291245	0.0	0.000183941874368	0.00187292768169	0.00164594001463	0.00162453700695" > exp.$$
gemini query --header -q "select chrom, start, end, ref, alt, aaf_gnomad_all, aaf_gnomad_all, \
	                        aaf_gnomad_afr, aaf_gnomad_amr, aaf_gnomad_eas, aaf_gnomad_fin, \
		                        aaf_gnomad_nfe, aaf_gnomad_oth, aaf_gnomad_sas from variants" test.exac.db > obs.$$


check obs.$$ exp.$$
#rm obs.$$ exp.$$

#############################################
