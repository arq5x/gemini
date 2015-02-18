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
# 1. Test set_somatic
####################################################################
echo -e "    somatic.t01...\c"
echo "H_LS-E2-A14P-01A-31D-A19H-09	C/T	0.1875	12	64	H_LS-E2-A14P-10A-01D-A19H-09	C/C	0.0	0	30	chr2	128046288	128046289	C	T	ERCC3
H_LS-E2-A14P-01A-31D-A19H-09	C/A	0.555555555556	10	18	H_LS-E2-A14P-10A-01D-A19H-09	C/C	0.0	0	17	chr17	7578460	7578461	C	A	TP53
H_LS-E2-A14P-01A-31D-A19H-09	CTGCTATTTTG/CG	0.22	11	50	H_LS-E2-A14P-10A-01D-A19H-09	CTGCTATTTTG/CTGCTATTTTG	0.0	0	70	chr17	59861630	59861641	CTGCTATTTTG	CG	BRIP1
H_LS-E2-A14P-01A-31D-A19H-09	GTTTTTTTTTTTTTGTATTTTTAGTAG/GTTTTTTTTTTTTGTATTTTTAGTAG	0.166666666667	4	24	H_LS-E2-A14P-10A-01D-A19H-09	GTTTTTTTTTTTTTGTATTTTTAGTAG/GTTTTTTTTTTTTTGTATTTTTAGTAG	0.0833333333333	1	12	chr6	132856479	132856506	GTTTTTTTTTTTTTGTATTTTTAGTAG	GTTTTTTTTTTTTGTATTTTTAGTAG	TAAR9
H_LS-E2-A14P-01A-31D-A19H-09	T/C	0.0972222222222	7	72	H_LS-E2-A14P-10A-01D-A19H-09	T/T	0.0697674418605	3	43	chr6	132857168	132857169	T	C	TAAR9
H_LS-E2-A14P-01A-31D-A19H-09	A/G	0.0969696969697	16	165	H_LS-E2-A14P-10A-01D-A19H-09	A/A	0.00653594771242	1	153	chr6	132922237	132922238	A	G	None
H_LS-E2-A14P-01A-31D-A19H-09	GAAAAAAAAAAAAAGGTGAAAATT/GAAAAAAAAAAAAGGTGAAAATT	0.217391304348	5	23	H_LS-E2-A14P-10A-01D-A19H-09	GAAAAAAAAAAAAAGGTGAAAATT/GAAAAAAAAAAAAAGGTGAAAATT	0.0	0	25	chrX	132838304	132838328	GAAAAAAAAAAAAAGGTGAAAATT	GAAAAAAAAAAAAGGTGAAAATT	GPC3
Identified and set 7 somatic mutations" > exp

gemini set_somatic test.somatic.db > obs

check obs exp
rm obs exp

####################################################################
# 2. Test actionable_mutations
####################################################################
echo -e "    somatic.t02...\c"

echo "tum_name	chrom	start	end	ref	alt	gene	impact	is_somatic	in_cosmic_census	dgidb_info
H_LS-E2-A14P-01A-31D-A19H-09	chr2	128046288	128046289	C	T	ERCC3	non_syn_coding	1	1	None
H_LS-E2-A14P-01A-31D-A19H-09	chr17	7578460	7578461	C	A	TP53	non_syn_coding	1	1	{'searchTerm': 'TP53', 'geneCategories': ['CLINICALLY ACTIONABLE', 'DRUGGABLE GENOME', 'TUMOR SUPPRESSOR', 'TRANSCRIPTION FACTOR COMPLEX', 'DRUG RESISTANCE', 'HISTONE MODIFICATION', 'DNA REPAIR', 'TRANSCRIPTION FACTOR BINDING'], 'geneName': 'TP53', 'geneLongName': 'tumor protein p53', 'interactions': [{'source': 'DrugBank', 'interactionId': '711cbe42-4930-4b46-963e-79ab35bbbd0f', 'interactionType': 'n/a', 'drugName': '1-(9-ETHYL-9H-CARBAZOL-3-YL)-N-METHYLMETHANAMINE'}, {'source': 'PharmGKB', 'interactionId': '8234d9b9-085d-49b1-aac2-cf5375d91477', 'interactionType': 'n/a', 'drugName': 'FLUOROURACIL'}, {'source': 'PharmGKB', 'interactionId': '605d7bca-7ed9-428e-aa7c-f76aafd66b54', 'interactionType': 'n/a', 'drugName': 'PACLITAXEL'}, {'source': 'TTD', 'interactionId': '1fe9db63-3581-435b-b22a-12d45c8c9864', 'interactionType': 'activator', 'drugName': 'CURAXIN CBLC102'}, {'source': 'TALC', 'interactionId': '8f8f6822-cb9e-40aa-8360-5532e059f1e7', 'interactionType': 'vaccine', 'drugName': 'EP-2101'}, {'source': 'TALC', 'interactionId': 'd59e14bc-b9a5-4c9f-a5aa-7ba322f0fa0e', 'interactionType': 'vaccine', 'drugName': 'MUTANT P53 PEPTIDE PULSED DENDRITIC CELL'}, {'source': 'TALC', 'interactionId': '79256b6e-9a16-4fbe-a237-28dbca28bc2a', 'interactionType': 'vaccine', 'drugName': 'AD.P53-DC'}]}
H_LS-E2-A14P-01A-31D-A19H-09	chr17	59861630	59861641	CTGCTATTTTG	CG	BRIP1	inframe_codon_loss	1	1	None
H_LS-E2-A14P-01A-31D-A19H-09	chrX	132838304	132838328	GAAAAAAAAAAAAAGGTGAAAATT	GAAAAAAAAAAAAGGTGAAAATT	GPC3	splice_region	1	1	None" > exp

gemini actionable_mutations test.somatic.db > obs

check obs exp
rm obs exp
