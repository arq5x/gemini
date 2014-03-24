#################################################################
# 1. Test pathway tool for VEP (75) annotations
#################################################################
echo "    path.t01...\c"
echo "chrom	start	end	ref	alt	impact	sample	genotype	gene	transcript	pathway
chr10	52004314	52004315	T	C	intron	M10500	C/C	ASAH2	ENST00000329428	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	52004314	52004315	T	C	intron	M128215	C/C	ASAH2	ENST00000329428	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	52004314	52004315	T	C	intron	M10500	C/C	ASAH2	ENST00000447815	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	52004314	52004315	T	C	intron	M128215	C/C	ASAH2	ENST00000447815	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	52004314	52004315	T	C	intron	M10500	C/C	ASAH2	ENST00000395526	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	52004314	52004315	T	C	intron	M128215	C/C	ASAH2	ENST00000395526	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000531469	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000309035	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000494626	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000337195	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000411419	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	135336655	135336656	G	A	intron	M10478	A/A	CYP2E1	ENST00000463117	hsa00982:Drug_metabolism_cytochrome_P450,hsa01100:Metabolic_pathways,hsa00590:Arachidonic_acid_metabolism,hsa00980:Metabolism_of_xenobiotics_by_cytochrome_P450,hsa00591:Linoleic_acid_metabolism
chr10	135336655	135336656	G	A	intron	M128215	A/A	CYP2E1	ENST00000463117	hsa00982:Drug_metabolism_cytochrome_P450,hsa01100:Metabolic_pathways,hsa00590:Arachidonic_acid_metabolism,hsa00980:Metabolism_of_xenobiotics_by_cytochrome_P450,hsa00591:Linoleic_acid_metabolism
chr10	135336655	135336656	G	A	upstream	M10478	A/A	CYP2E1	ENST00000252945	hsa00982:Drug_metabolism_cytochrome_P450,hsa01100:Metabolic_pathways,hsa00590:Arachidonic_acid_metabolism,hsa00980:Metabolism_of_xenobiotics_by_cytochrome_P450,hsa00591:Linoleic_acid_metabolism
chr10	135336655	135336656	G	A	upstream	M128215	A/A	CYP2E1	ENST00000252945	hsa00982:Drug_metabolism_cytochrome_P450,hsa01100:Metabolic_pathways,hsa00590:Arachidonic_acid_metabolism,hsa00980:Metabolism_of_xenobiotics_by_cytochrome_P450,hsa00591:Linoleic_acid_metabolism
chr16	72057434	72057435	C	T	non_syn_coding	M10475	C/T	DHODH	ENST00000219240	hsa01100:Metabolic_pathways,hsa00240:Pyrimidine_metabolism" > exp
gemini pathways -v 68 test4.vep.db \
       > obs
check obs exp
rm obs exp

#################################################################
# 2. Test pathway tool with -lof option for VEP annotations
#################################################################
echo "    path.t02...\c"
echo "chrom	start	end	ref	alt	impact	sample	genotype	gene	transcript	pathway
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000531469	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000309035	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000494626	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000337195	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000411419	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer" > exp
gemini pathways -v 68 --lof test4.vep.db \
       > obs
check obs exp
rm obs exp

#################################################################
# 3. Test pathway tool for snpEff (v3.0) annotations
#################################################################
echo "    path.t03...\c"
echo "chrom	start	end	ref	alt	impact	sample	genotype	gene	transcript	pathway
chr10	52004314	52004315	T	C	intron	M10500	C/C	ASAH2	ENST00000395526	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	52004314	52004315	T	C	intron	M128215	C/C	ASAH2	ENST00000395526	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000337195	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000411419	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	135336655	135336656	G	A	upstream	M10478	A/A	CYP2E1	ENST00000252945	hsa00982:Drug_metabolism_cytochrome_P450,hsa01100:Metabolic_pathways,hsa00590:Arachidonic_acid_metabolism,hsa00980:Metabolism_of_xenobiotics_by_cytochrome_P450,hsa00591:Linoleic_acid_metabolism
chr10	135336655	135336656	G	A	upstream	M128215	A/A	CYP2E1	ENST00000252945	hsa00982:Drug_metabolism_cytochrome_P450,hsa01100:Metabolic_pathways,hsa00590:Arachidonic_acid_metabolism,hsa00980:Metabolism_of_xenobiotics_by_cytochrome_P450,hsa00591:Linoleic_acid_metabolism
chr16	72057434	72057435	C	T	non_syn_coding	M10475	C/T	DHODH	ENST00000219240	hsa01100:Metabolic_pathways,hsa00240:Pyrimidine_metabolism" > exp
gemini pathways -v 66 test4.snpeff.db \
       > obs
check obs exp
rm obs exp

#################################################################
# 4. Test pathway tool for --lof option with snpEff (v3.0) annotations
#################################################################
echo "    path.t04...\c"
echo "chrom	start	end	ref	alt	impact	sample	genotype	gene	transcript	pathway
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000337195	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000411419	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer" > exp
gemini pathways -v 66 --lof test4.snpeff.db \
       > obs
check obs exp
rm obs exp

