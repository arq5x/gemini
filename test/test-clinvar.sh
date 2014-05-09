
####################################################################
# 1. Test clinvar annotations
####################################################################
echo "    clinvar.t01...\c"
echo "1	pathogenic	Myasthenia,_limb-girdle,_familial	OMIM_Allelic_Variant	103320.0001	germline	GeneReviews:MedGen:OMIM:Orphanet	NBK1168:C1850792:254300:590	RCV000019902.26	1	0
None	None	None	None	None	None	None	None	None	None	None
1	other	Generalized_epilepsy_with_febrile_seizures_plus_type_5|Epilepsy,_juvenile_myoclonic_7|Epilepsy,_idiopathic_generalized_10	OMIM_Allelic_Variant	137163.0002	germline	MedGen|MedGen|MedGen:OMIM	C3150401|CN043549|C2751603:613060	RCV000017599.1|RCV000017600.1|RCV000022558.1	1	0
1	pathogenic	Roussy-Levy_syndrome	OMIM_Allelic_Variant	159440.0021	germline	MedGen:OMIM:SNOMED_CT	C0205713:180800:45853006	RCV000015250.24	1	0
1	pathogenic	Chediak-Higashi_syndrome	GeneReviews	NBK5188	unknown	GeneReviews:MedGen:OMIM:Orphanet:SNOMED_CT	NBK5188:C0007965:214500:167:111396008	RCV000033871.2	1	0
1	untested	Familial_cold_urticaria	Unite_medicale_des_maladies_autoinflammatoires	363	None	MedGen:OMIM:Orphanet:SNOMED_CT	C0343068:120100:47045:238687000	RCV000084222.1	1	0
1	untested	Juvenile_GM>1<_gangliosidosis	.	.	somatic	MedGen:OMIM:Orphanet:Orphanet:SNOMED_CT	C0268272:230600:354:79256:18756002	RCV000056404.1	1	0" > exp

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