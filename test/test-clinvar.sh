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
echo "chr1	985955	1	pathogenic	Myasthenic_syndrome,_congenital,_with_pre-_and_postsynaptic_defects	OMIM_Allelic_Variant	103320.0001	germline	MedGen:OMIM:Orphanet	CN186031:615120:ORPHA590	RCV000019902.27	1	0	C
chr1	1199489	1	other	Lung_cancer	.	.	somatic	MedGen:OMIM:SNOMED_CT	C0684249:211980:187875007	RCV000089747.1	1	0	A
chr1	1959699	1	other	Generalized_epilepsy_with_febrile_seizures_plus_type_5|Epilepsy,_juvenile_myoclonic_7|Epilepsy,_idiopathic_generalized_10	OMIM_Allelic_Variant	137163.0002	germline	MedGen|MedGen|MedGen:OMIM	C3150401|CN043549|C2751603:613060	RCV000017599.1|RCV000017600.1|RCV000022558.1	1	0	A
chr1	161276553	1	pathogenic	Roussy-Levy_syndrome	OMIM_Allelic_Variant	159440.0021	germline	MedGen:OMIM:SNOMED_CT	C0205713:180800:45853006	RCV000015250.24	1	0	T
chr1	235891431	None	None	None	None	None	None	None	None	None	None	None	None
chr1	247587093	1	not-provided	Familial_cold_urticaria	.	.	unknown	MedGen:OMIM:Orphanet:SNOMED_CT	C0343068:120100:ORPHA47045:238687000	RCV000084222.1	1	0	T
chr3	33063141	1	other	Juvenile_GM>1<_gangliosidosis	.	.	somatic	MedGen:OMIM:Orphanet:Orphanet:SNOMED_CT	C0268272:230600:ORPHA354:ORPHA79256:18756002	RCV000056404.1	1	0	A
chrX	31200832	0	benign	not_specified,not_specified	.	.	None	MedGen,MedGen	CN169374,CN169374	RCV000080866.2,RCV000080867.2	1	0	A,AT
chrX	89963314	None	None	None	None	None	None	None	None	None	None	None	None" > exp

gemini query -q "select chrom, end, in_omim,
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
