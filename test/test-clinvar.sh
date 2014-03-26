
####################################################################
# 1. Test clinvar annotations
####################################################################
echo "    clinvar.t01...\c"
echo "1	pathogenic	Myasthenia,_limb-girdle,_familial	OMIM_Allelic_Variant	103320.0001	germline	GeneReviews:MedGen:OMIM:Orphanet	NBK1168:C1850792:254300:590	RCV000019902.26	1	0
None	None	None	None	None	None	None	None	None	None	None
0	untested	Malignant_melanoma	.	.	somatic	MedGen:SNOMED_CT	C0025202:2092003	RCV000059856.2	1	0
0	untested	Malignant_melanoma	.	.	somatic	MedGen:SNOMED_CT	C0025202:2092003	RCV000059865.2	1	0
0	untested	Malignant_melanoma	.	.	somatic	MedGen:SNOMED_CT	C0025202:2092003	RCV000064284.2	1	0
None	None	None	None	None	None	None	None	None	None	None
1	other	Generalized_epilepsy_with_febrile_seizures_plus_type_5	OMIM_Allelic_Variant	137163.0001	germline	MedGen	C3150401	RCV000017598.1	1	0
1	other	Generalized_epilepsy_with_febrile_seizures_plus_type_5|Epilepsy,_juvenile_myoclonic_7|Epilepsy,_idiopathic_generalized_10	OMIM_Allelic_Variant	137163.0002	germline	MedGen|MedGen|MedGen:OMIM	C3150401|CN043549|C2751603:613060	RCV000017599.1|RCV000017600.1|RCV000022558.1	1	0
0	untested	Malignant_melanoma	.	.	somatic	MedGen:SNOMED_CT	C0025202:2092003	RCV000060029.2	1	0
None	None	None	None	None	None	None	None	None	None	None" > exp

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