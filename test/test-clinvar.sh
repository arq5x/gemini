check()
{
        if diff $1 $2; then
        echo ok
        else
        echo fail
        fi
}
export -f check
####################################################################
# 1. Test clinvar annotations
####################################################################
echo "    clinvar.t01...\c"
echo "1	pathogenic	MYASTHENIC_SYNDROME,_CONGENITAL,_WITH_PRE-_AND_POSTSYNAPTIC_DEFECTS	ClinVar|OMIM_Allelic_Variant	NM_198576.3:c.5125G>C|103320.0001	germline	MedGen:OMIM:Orphanet	CN186031:615120:ORPHA590	RCV000019902.27	1	0	C
1	other	Lung_cancer	ClinVar	NM_058167.2:c.132-723C>T	somatic	GeneReviews:MedGen:OMIM:SNOMED_CT	NBK1108:C0684249:211980:187875007	RCV000089747.1	1	0	A
1	other	Generalized_epilepsy_with_febrile_seizures_plus_type_5|Epilepsy,_juvenile_myoclonic_7|Epilepsy,_idiopathic_generalized_10	ClinVar|OMIM_Allelic_Variant	NM_000815.4:c.659G>A|137163.0002	germline	MedGen|MedGen|MedGen:OMIM	C3150401|CN043549|C2751603:613060	RCV000017599.1|RCV000017600.1|RCV000022558.1	1	0	A
1	pathogenic	Roussy-Levy_syndrome	ClinVar|OMIM_Allelic_Variant	NM_000530.6:c.393C>A|159440.0021	germline	MedGen:OMIM:SNOMED_CT	C0205713:180800:45853006	RCV000015250.24	1	0	T
1	pathogenic	Chediak-Higashi_syndrome	ClinVar|GeneReviews	NM_000081.3:c.9107_9162del56|NBK5188	unknown	GeneReviews:MedGen:OMIM:Orphanet:SNOMED_CT	NBK5188:C0007965:214500:ORPHA167:111396008	RCV000033871.2	1	0	C
1	untested	Familial_cold_urticaria	ClinVar|Unite_medicale_des_maladies_autoinflammatoires	NM_004895.4:c.404-56C>T|363	unknown	MedGen:OMIM:Orphanet:SNOMED_CT	C0343068:120100:ORPHA47045:238687000	RCV000084222.1	1	0	T
1	other	Juvenile_GM>1<_gangliosidosis	ClinVar	NM_000404.2:c.1150A>T	somatic	MedGen:OMIM:Orphanet:Orphanet:SNOMED_CT	C0268272:230600:ORPHA354:ORPHA79256:18756002	RCV000056404.1	1	0	A
0	None	AllHighlyPenetrant,AllHighlyPenetrant	.,.	.,.	None	.,.	.,.	RCV000080866.1,RCV000080867.1	1	0	A,AT
1	other	Lung_cancer	.	.	somatic	GeneReviews:MedGen:OMIM:SNOMED_CT	NBK1108:C0684249:211980:187875007	RCV000104117.1	1	0	T" > exp

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
                        clinvar_on_diag_assay,
                        clinvar_causal_allele from variants" test.clinvar.db > obs
check obs exp
#rm obs exp