
####################################################################
# 1. 
####################################################################
echo "    clinvar.t01...\c"
echo "1	pathogenic	Myasthenia, limb-girdle, familial	OMIM Allelic Variant	103320.0001	germline	GeneReviews:NCBI:OMIM:Orphanet	NBK1168:C1850792:254300:590	RCV000019902.1	1	0
1	untested	.	.	.	somatic	.	.	.	1	0
0	untested	.	.	.	somatic	.	.	.	1	0
0	untested	.	.	.	somatic	.	.	.	1	0
0	untested	.	.	.	somatic	.	.	.	1	0
1	untested	.	.	.	somatic	.	.	.	1	0
1	other	Generalized epilepsy with febrile seizures plus type 5	OMIM Allelic Variant	137163.0001	germline	NCBI	C1858675	RCV000017598.1	1	0
1	other	Generalized epilepsy with febrile seizures plus type 5|Epilepsy, juvenile myoclonic 7|Epilepsy, idiopathic generalized 10	OMIM Allelic Variant	137163.0002	germline	NCBI|NCBI|NCBI:OMIM	C1858675|CN043549|C2751603:613060	RCV000017599.1|RCV000017600.1|RCV000022558.1	1	0
0	untested	.	.	.	somatic	.	.	.	1	0
0	untested	.	.	.	unknown	.	.	.	1	0" > exp

gemini query -q "select in_omim, 
                        clinvar_sig, 
                        clinvar_disease_name, 
                        clinvar_dbsource, 
                        clinvar_dbsource_id, 
                        clinvar_origin, 
                        clinvar_dsdb, 
                        clinvar_dsdbid, 
                        clinvar_disease_acc, 
                        clinvar_in_locus_spec_db, 
                        clinvar_on_diag_assay from variants" test.clinvar.db > obs
check obs exp
rm obs exp
