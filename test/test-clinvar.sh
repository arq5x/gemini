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
echo "    clinvar.t01..."
echo "chr1	985955	1	pathogenic	Congenital_myasthenic_syndrome|Myasthenic_syndrome,_congenital,_8	germline	MedGen:C0751882,Orphanet:ORPHA590|MedGen:C3808739,OMIM:615120	None
chr1	1199489	0	None	None	None	None	None
chr1	1959699	1	benign/likely_benign,_risk_factor	Idiopathic_generalized_epilepsy|Epilepsy,_idiopathic_generalized_10|Epilepsy,_juvenile_myoclonic_7|Generalized_epilepsy_with_febrile_seizures_plus_type_5|not_provided	germline	Gene:1957,MedGen:C0270850,OMIM:600669,SNOMED_CT:36803009|MedGen:C2751603,OMIM:613060|MedGen:C2751604|MedGen:C3150401|MedGen:CN517202	None
chr1	161276553	1	pathogenic	Roussy-Levy_syndrome|Charcot-Marie-Tooth_disease,_demyelinating,_type_1b|not_provided	germline	MedGen:C0205713,OMIM:180800,Orphanet:ORPHA3115,SNOMED_CT:45853006|MedGen:C0270912,OMIM:118200,Orphanet:ORPHA101082,SNOMED_CT:42986003|MedGen:CN517202	None
chr1	235891431	0	None	None	None	None	None
chr1	247587093	1	not_provided	Familial_cold_urticaria	unknown	MedGen:C0343068,OMIM:120100,Orphanet:ORPHA47045,SNOMED_CT:238687000	None
chr3	33063141	1	not_provided	GM1_gangliosidosis_type_2	somatic	MedGen:C0268272,OMIM:230600,Orphanet:ORPHA79256,SNOMED_CT:18756002	None
chrX	31200832	0	None	None	None	None	None
chrX	89963314	0	None	None	None	None	None" > exp
gemini query -q "select chrom, end, in_omim,
                        clinvar_sig, 
                        clinvar_disease_name, 
                        clinvar_origin, 
                        clinvar_dsdb, 
                        clinvar_dsdbid from variants" test.clinvar.db > obs
check obs exp
rm obs exp

echo "    clinvar.t02..."
echo "B3GALT6	ehlers-danlos_syndrome,_progeroid_type,_2|spondyloepimetaphyseal_dysplasia_with_joint_laxity|spondyloepimetaphyseal_dysplasia_with_joint_laxity,_type_1,_with_fractures
SCNN1D	None" > exp

 gemini query -q "select gene, clinvar_gene_phenotype from variants" test.clinvar_gene_pheno.db > obs
 check obs exp
