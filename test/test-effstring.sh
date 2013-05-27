####################################################################
# 1. Test variant_impacts table for the impact columns
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
echo "anno_id	gene	transcript	impact	impact_severity	biotype	is_exonic	is_coding	is_lof
1	FAM138A	ENST00000417324	downstream	LOW	protein_coding	0	0	0
1	OR4F5	ENST00000335137	synonymous_coding	LOW	protein_coding	1	1	0
1	OR4F5	ENST00000335137	non_syn_coding	MED	protein_coding	1	1	0
6	SAMD11	ENST00000342066	frame_shift	HIGH	protein_coding	1	1	1
8	SAMD11	ENST00000342066	UTR_3_prime	LOW	protein_coding	1	0	0
1	HES4	ENST00000428771	start_gain	LOW	protein_coding	1	0	0
13	CCNL2	ENST00000488340	splice_acceptor	HIGH	nonsense_mediated_decay	0	0	0" > exp

gemini query -q "select anno_id, gene, transcript, impact, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variants" \
                    --header \
                    test1.snpeff.db \
                    > obs
check obs exp
rm obs exp
