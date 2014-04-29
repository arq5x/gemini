####################################################################
# 1. Test variant_impacts table for the impact columns (snpEff)
####################################################################
echo "    effstring.t01...\c"
echo "gene	impact	impact_severity	biotype	is_exonic	is_coding	is_lof
FAM138A	downstream	LOW	protein_coding	0	0	0
FAM138A	downstream	LOW	processed_transcript	0	0	0
MIR1302-10	downstream	LOW	miRNA	0	0	0
MIR1302-10	intron	LOW	antisense	0	0	0
MIR1302-10	intron	LOW	antisense	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
OR4F5	synonymous_coding	LOW	protein_coding	1	1	0
OR4F5	non_syn_coding	MED	protein_coding	1	1	0
SAMD11	downstream	LOW	protein_coding	0	0	0
SAMD11	downstream	LOW	protein_coding	0	0	0
NOC2L	downstream	LOW	protein_coding	0	0	0
NOC2L	downstream	LOW	processed_transcript	0	0	0
NOC2L	downstream	LOW	processed_transcript	0	0	0
SAMD11	frame_shift	HIGH	protein_coding	1	1	1
SAMD11	upstream	LOW	processed_transcript	0	0	0
SAMD11	upstream	LOW	processed_transcript	0	0	0
SAMD11	upstream	LOW	retained_intron	0	0	0
SAMD11	upstream	LOW	retained_intron	0	0	0
NOC2L	downstream	LOW	processed_transcript	0	0	0
SAMD11	downstream	LOW	processed_transcript	0	0	0
SAMD11	downstream	LOW	processed_transcript	0	0	0
SAMD11	downstream	LOW	retained_intron	0	0	0
SAMD11	downstream	LOW	retained_intron	0	0	0
NOC2L	exon	LOW	processed_transcript	0	0	0
NOC2L	exon	LOW	processed_transcript	0	0	0
SAMD11	UTR_3_prime	LOW	protein_coding	1	0	0
NOC2L	UTR_3_prime	LOW	protein_coding	1	0	0
HES4	start_gain	LOW	protein_coding	1	0	0
HES4	upstream	LOW	protein_coding	0	0	0
HES4	upstream	LOW	protein_coding	0	0	0
HES4	upstream	LOW	processed_transcript	0	0	0
MRPL20	downstream	LOW	protein_coding	0	0	0
MRPL20	downstream	LOW	processed_transcript	0	0	0
MRPL20	downstream	LOW	processed_transcript	0	0	0
MRPL20	downstream	LOW	processed_transcript	0	0	0
RP4-758J18.5.1	downstream	LOW	processed_transcript	0	0	0
CCNL2	exon	LOW	processed_transcript	0	0	0
CCNL2	intron	LOW	protein_coding	0	0	0
CCNL2	intron	LOW	protein_coding	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	splice_acceptor	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	upstream	LOW	processed_transcript	0	0	0
CCNL2	upstream	LOW	processed_transcript	0	0	0
CCNL2	upstream	LOW	processed_transcript	0	0	0
CCNL2	upstream	LOW	retained_intron	0	0	0
RP4-758J18.2.1	upstream	LOW	processed_transcript	0	0	0
RP4-758J18.2.1	upstream	LOW	processed_transcript	0	0	0
RP4-758J18.2.1	upstream	LOW	processed_transcript	0	0	0
RP4-758J18.2.1	upstream	LOW	processed_transcript	0	0	0
RP4-758J18.3.1	upstream	LOW	processed_transcript	0	0	0" > exp

gemini query -q "select gene, impact, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variant_impacts" \
                    --header \
                    test1.snpeff.db \
                    > obs
check obs exp
rm obs exp

####################################################################
# 2. Test variants table for severe_impact columns  
####################################################################
echo "    effstring.t02...\c"
echo "anno_id	gene	transcript	impact	impact_so	impact_severity	biotype	is_exonic	is_coding	is_lof
1	FAM138A	ENST00000417324	downstream	downstream_gene_variant	LOW	protein_coding	0	0	0
1	OR4F5	ENST00000335137	synonymous_coding	synonymous_variant	LOW	protein_coding	1	1	0
1	OR4F5	ENST00000335137	non_syn_coding	missense_variant	MED	protein_coding	1	1	0
6	SAMD11	ENST00000342066	frame_shift	frameshift_variant	HIGH	protein_coding	1	1	1
8	SAMD11	ENST00000342066	UTR_3_prime	3_prime_UTR_variant	LOW	protein_coding	1	0	0
1	HES4	ENST00000428771	start_gain	5_prime_UTR_premature_start_codon_gain_variant	LOW	protein_coding	1	0	0
13	CCNL2	ENST00000488340	splice_acceptor	splice_acceptor_variant	HIGH	nonsense_mediated_decay	0	0	0" > exp

gemini query -q "select anno_id, gene, transcript, impact, impact_so, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variants" \
                    --header \
                    test1.snpeff.db \
                    > obs
check obs exp
rm obs exp

###################################################################
# 3. Test variant_impacts table for the impact columns (VEP)
###################################################################
echo "    effstring.t03...\c"
echo "gene	impact	impact_severity	biotype	is_exonic	is_coding	is_lof
MIR1302-10	intron	LOW	lincRNA	0	0	0
MIR1302-10	nc_transcript	LOW	lincRNA	0	0	0
MIR1302-10	intron	LOW	lincRNA	0	0	0
MIR1302-10	nc_transcript	LOW	lincRNA	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
FAM138A	downstream	LOW	lincRNA	0	0	0
MIR1302-10	downstream	LOW	miRNA	0	0	0
FAM138A	downstream	LOW	lincRNA	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream	LOW	unprocessed_pseudogene	0	0	0
OR4F5	synonymous_coding	LOW	protein_coding	1	1	0
OR4F5	non_syn_coding	MED	protein_coding	1	1	0
SAMD11	downstream	LOW	protein_coding	0	0	0
SAMD11	frame_shift	HIGH	protein_coding	1	1	1
SAMD11	feature_elongation	LOW	protein_coding	1	1	0
NOC2L	downstream	LOW	protein_coding	0	0	0
SAMD11	upstream	LOW	retained_intron	0	0	0
SAMD11	upstream	LOW	retained_intron	0	0	0
SAMD11	upstream	LOW	retained_intron	0	0	0
NOC2L	downstream	LOW	retained_intron	0	0	0
SAMD11	frame_shift	HIGH	protein_coding	1	1	1
SAMD11	feature_elongation	LOW	protein_coding	1	1	0
SAMD11	downstream	LOW	protein_coding	0	0	0
SAMD11	upstream	LOW	processed_transcript	0	0	0
NOC2L	downstream	LOW	retained_intron	0	0	0
SAMD11	intron	LOW	protein_coding	0	0	0
SAMD11	feature_elongation	LOW	protein_coding	1	1	0
NOC2L	downstream	LOW	processed_transcript	0	0	0
SAMD11	UTR_3_prime	LOW	protein_coding	1	0	0
NOC2L	UTR_3_prime	LOW	protein_coding	1	0	0
SAMD11	downstream	LOW	retained_intron	0	0	0
SAMD11	downstream	LOW	retained_intron	0	0	0
SAMD11	downstream	LOW	retained_intron	0	0	0
SAMD11	downstream	LOW	protein_coding	0	0	0
NOC2L	nc_exon	LOW	retained_intron	0	0	0
NOC2L	nc_transcript	LOW	retained_intron	0	0	0
SAMD11	downstream	LOW	processed_transcript	0	0	0
NOC2L	nc_exon	LOW	retained_intron	0	0	0
NOC2L	nc_transcript	LOW	retained_intron	0	0	0
SAMD11	UTR_3_prime	LOW	protein_coding	1	0	0
HES4	upstream	LOW	protein_coding	0	0	0
HES4	UTR_5_prime	LOW	protein_coding	1	0	0
RP11-54O7.17	upstream	LOW	lincRNA	0	0	0
HES4	upstream	LOW	retained_intron	0	0	0
HES4	upstream	LOW	protein_coding	0	0	0
RP4-758J18.2	upstream	LOW	protein_coding	0	0	0
CCNL2	splice_acceptor	HIGH	protein_coding	0	0	1
CCNL2	intron	LOW	protein_coding	0	0	0
CCNL2	feature_truncation	LOW	protein_coding	1	1	0
RP4-758J18.2	upstream	LOW	retained_intron	0	0	0
RP4-758J18.2	upstream	LOW	protein_coding	0	0	0
CCNL2	upstream	LOW	retained_intron	0	0	0
CCNL2	upstream	LOW	processed_transcript	0	0	0
CCNL2	splice_acceptor	HIGH	protein_coding	0	0	1
CCNL2	intron	LOW	protein_coding	0	0	0
CCNL2	feature_truncation	LOW	protein_coding	1	1	0
MRPL20	downstream	LOW	retained_intron	0	0	0
RP4-758J18.2	upstream	LOW	processed_transcript	0	0	0
MRPL20	downstream	LOW	protein_coding	0	0	0
RP4-758J18.2	upstream	LOW	retained_intron	0	0	0
CCNL2	upstream	LOW	retained_intron	0	0	0
CCNL2	splice_acceptor	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	NMD_transcript	LOW	nonsense_mediated_decay	0	0	0
CCNL2	feature_truncation	LOW	nonsense_mediated_decay	0	0	0
CCNL2	upstream	LOW	nonsense_mediated_decay	0	0	0
RP4-758J18.2	upstream	LOW	protein_coding	0	0	0
CCNL2	splice_acceptor	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	NMD_transcript	LOW	nonsense_mediated_decay	0	0	0
CCNL2	feature_truncation	LOW	nonsense_mediated_decay	0	0	0
CCNL2	UTR_5_prime	LOW	protein_coding	1	0	0
CCNL2	feature_truncation	LOW	protein_coding	1	1	0
CCNL2	upstream	LOW	protein_coding	0	0	0
CCNL2	splice_acceptor	HIGH	retained_intron	0	0	0
CCNL2	intron	LOW	retained_intron	0	0	0
CCNL2	nc_transcript	LOW	retained_intron	0	0	0
CCNL2	feature_truncation	LOW	retained_intron	0	0	0
RP4-758J18.2	upstream	LOW	protein_coding	0	0	0
CCNL2	splice_acceptor	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	intron	LOW	nonsense_mediated_decay	0	0	0
CCNL2	NMD_transcript	LOW	nonsense_mediated_decay	0	0	0
CCNL2	feature_truncation	LOW	nonsense_mediated_decay	0	0	0
MRPL20	downstream	LOW	processed_transcript	0	0	0
RP4-758J18.2	upstream	LOW	nonsense_mediated_decay	0	0	0
MRPL20	downstream	LOW	protein_coding	0	0	0" > exp

gemini query -q "select gene, impact, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variant_impacts" \
                    --header \
                    test1.vep.db \
                    > obs
check obs exp
rm obs exp

#########################################################################
# 4. Test variants table for severe_impact columns (VEP)
#########################################################################
echo "    effstring.t04...\c"
echo "anno_id	gene	transcript	impact	impact_severity	biotype	is_exonic	is_coding	is_lof
1	MIR1302-10	ENST00000473358	intron	LOW	lincRNA	0	0	0
1	OR4F5	ENST00000335137	synonymous_coding	LOW	protein_coding	1	1	0
1	OR4F5	ENST00000335137	non_syn_coding	MED	protein_coding	1	1	0
2	SAMD11	ENST00000342066	frame_shift	HIGH	protein_coding	1	1	1
2	SAMD11	ENST00000342066	UTR_3_prime	LOW	protein_coding	1	0	0
1	HES4	ENST00000484667	upstream	LOW	protein_coding	0	0	0
2	CCNL2	ENST00000400809	splice_acceptor	HIGH	protein_coding	0	0	1" > exp

gemini query -q "select anno_id, gene, transcript, impact, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variants" \
                    --header \
                    test1.vep.db \
                    > obs
check obs exp
rm obs exp

###########################################################################