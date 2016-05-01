####################################################################
# 1. Test variant_impacts table for the impact columns (snpEff)
####################################################################
echo "    effstring.t01...\c"
echo "gene	impact	impact_severity	biotype	is_exonic	is_coding	is_lof
FAM138A	downstream_gene_variant	LOW	protein_coding	0	0	0
FAM138A	downstream_gene_variant	LOW	processed_transcript	0	0	0
MIR1302-10	downstream_gene_variant	LOW	miRNA	0	0	0
MIR1302-10	intron_variant	LOW	antisense	0	0	0
MIR1302-10	intron_variant	LOW	antisense	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
OR4F5	synonymous_variant	LOW	protein_coding	1	1	0
OR4F5	missense_variant	MED	protein_coding	1	1	0
SAMD11	downstream_gene_variant	LOW	protein_coding	0	0	0
SAMD11	downstream_gene_variant	LOW	protein_coding	0	0	0
NOC2L	downstream_gene_variant	LOW	protein_coding	0	0	0
NOC2L	downstream_gene_variant	LOW	processed_transcript	0	0	0
NOC2L	downstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	frameshift_variant	HIGH	protein_coding	1	1	1
SAMD11	upstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	upstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	upstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	upstream_gene_variant	LOW	retained_intron	0	0	0
NOC2L	downstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	downstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	downstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	downstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	downstream_gene_variant	LOW	retained_intron	0	0	0
NOC2L	exon_variant	LOW	processed_transcript	1	0	0
NOC2L	exon_variant	LOW	processed_transcript	1	0	0
SAMD11	3_prime_UTR_variant	LOW	protein_coding	1	0	0
NOC2L	3_prime_UTR_variant	LOW	protein_coding	1	0	0
HES4	5_prime_UTR_premature_start_codon_variant	LOW	protein_coding	1	0	0
HES4	upstream_gene_variant	LOW	protein_coding	0	0	0
HES4	upstream_gene_variant	LOW	protein_coding	0	0	0
HES4	upstream_gene_variant	LOW	processed_transcript	0	0	0
MRPL20	downstream_gene_variant	LOW	protein_coding	0	0	0
MRPL20	downstream_gene_variant	LOW	processed_transcript	0	0	0
MRPL20	downstream_gene_variant	LOW	processed_transcript	0	0	0
MRPL20	downstream_gene_variant	LOW	processed_transcript	0	0	0
RP4-758J18.5.1	downstream_gene_variant	LOW	processed_transcript	0	0	0
CCNL2	exon_variant	LOW	processed_transcript	1	0	0
CCNL2	intron_variant	LOW	protein_coding	0	0	0
CCNL2	intron_variant	LOW	protein_coding	0	0	0
CCNL2	intron_variant	LOW	nonsense_mediated_decay	0	0	0
CCNL2	intron_variant	LOW	nonsense_mediated_decay	0	0	0
CCNL2	intron_variant	LOW	nonsense_mediated_decay	0	0	0
CCNL2	intron_variant	LOW	nonsense_mediated_decay	0	0	0
CCNL2	splice_acceptor_variant	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	upstream_gene_variant	LOW	processed_transcript	0	0	0
CCNL2	upstream_gene_variant	LOW	processed_transcript	0	0	0
CCNL2	upstream_gene_variant	LOW	processed_transcript	0	0	0
CCNL2	upstream_gene_variant	LOW	retained_intron	0	0	0
RP4-758J18.2.1	upstream_gene_variant	LOW	processed_transcript	0	0	0
RP4-758J18.2.1	upstream_gene_variant	LOW	processed_transcript	0	0	0
RP4-758J18.2.1	upstream_gene_variant	LOW	processed_transcript	0	0	0
RP4-758J18.2.1	upstream_gene_variant	LOW	processed_transcript	0	0	0
RP4-758J18.3.1	upstream_gene_variant	LOW	processed_transcript	0	0	0" > exp

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
1	FAM138A	ENST00000417324	downstream_gene_variant	downstream_gene_variant	LOW	protein_coding	0	0	0
1	OR4F5	ENST00000335137	synonymous_variant	synonymous_variant	LOW	protein_coding	1	1	0
1	OR4F5	ENST00000335137	missense_variant	missense_variant	MED	protein_coding	1	1	0
6	SAMD11	ENST00000342066	frameshift_variant	frameshift_variant	HIGH	protein_coding	1	1	1
8	SAMD11	ENST00000342066	3_prime_UTR_variant	3_prime_UTR_variant	LOW	protein_coding	1	0	0
1	HES4	ENST00000428771	5_prime_UTR_premature_start_codon_variant	5_prime_UTR_premature_start_codon_variant	LOW	protein_coding	1	0	0
13	CCNL2	ENST00000488340	splice_acceptor_variant	splice_acceptor_variant	HIGH	nonsense_mediated_decay	0	0	0" > exp

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
echo "CCNL2	5_prime_UTR_variant	LOW	protein_coding	1	0	0
CCNL2	splice_acceptor_variant	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	splice_acceptor_variant	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	splice_acceptor_variant	HIGH	nonsense_mediated_decay	0	0	0
CCNL2	splice_acceptor_variant	HIGH	protein_coding	0	0	1
CCNL2	splice_acceptor_variant	HIGH	protein_coding	0	0	1
CCNL2	splice_acceptor_variant	HIGH	retained_intron	0	0	0
CCNL2	upstream_gene_variant	LOW	nonsense_mediated_decay	0	0	0
CCNL2	upstream_gene_variant	LOW	processed_transcript	0	0	0
CCNL2	upstream_gene_variant	LOW	protein_coding	0	0	0
CCNL2	upstream_gene_variant	LOW	retained_intron	0	0	0
CCNL2	upstream_gene_variant	LOW	retained_intron	0	0	0
FAM138A	downstream_gene_variant	LOW	lincRNA	0	0	0
FAM138A	downstream_gene_variant	LOW	lincRNA	0	0	0
HES4	5_prime_UTR_variant	LOW	protein_coding	1	0	0
HES4	upstream_gene_variant	LOW	protein_coding	0	0	0
HES4	upstream_gene_variant	LOW	protein_coding	0	0	0
HES4	upstream_gene_variant	LOW	retained_intron	0	0	0
MIR1302-10	downstream_gene_variant	LOW	miRNA	0	0	0
MIR1302-10	intron_variant	LOW	lincRNA	0	0	0
MIR1302-10	intron_variant	LOW	lincRNA	0	0	0
MRPL20	downstream_gene_variant	LOW	processed_transcript	0	0	0
MRPL20	downstream_gene_variant	LOW	protein_coding	0	0	0
MRPL20	downstream_gene_variant	LOW	protein_coding	0	0	0
MRPL20	downstream_gene_variant	LOW	retained_intron	0	0	0
NOC2L	3_prime_UTR_variant	LOW	protein_coding	1	0	0
NOC2L	downstream_gene_variant	LOW	processed_transcript	0	0	0
NOC2L	downstream_gene_variant	LOW	protein_coding	0	0	0
NOC2L	downstream_gene_variant	LOW	retained_intron	0	0	0
NOC2L	downstream_gene_variant	LOW	retained_intron	0	0	0
NOC2L	non_coding_exon_variant	LOW	retained_intron	0	0	0
NOC2L	non_coding_exon_variant	LOW	retained_intron	0	0	0
OR4F5	missense_variant	MED	protein_coding	1	1	0
OR4F5	synonymous_variant	LOW	protein_coding	1	1	0
RP11-54O7.17	upstream_gene_variant	LOW	lincRNA	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	nonsense_mediated_decay	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	processed_transcript	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	protein_coding	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	protein_coding	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	protein_coding	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	protein_coding	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	retained_intron	0	0	0
RP4-758J18.2	upstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	3_prime_UTR_variant	LOW	protein_coding	1	0	0
SAMD11	3_prime_UTR_variant	LOW	protein_coding	1	0	0
SAMD11	downstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	downstream_gene_variant	LOW	protein_coding	0	0	0
SAMD11	downstream_gene_variant	LOW	protein_coding	0	0	0
SAMD11	downstream_gene_variant	LOW	protein_coding	0	0	0
SAMD11	downstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	downstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	downstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	frameshift_variant	HIGH	protein_coding	1	1	1
SAMD11	frameshift_variant	HIGH	protein_coding	1	1	1
SAMD11	intron_variant	LOW	protein_coding	0	0	0
SAMD11	upstream_gene_variant	LOW	processed_transcript	0	0	0
SAMD11	upstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	upstream_gene_variant	LOW	retained_intron	0	0	0
SAMD11	upstream_gene_variant	LOW	retained_intron	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
WASH7P	upstream_gene_variant	LOW	unprocessed_pseudogene	0	0	0
gene	impact	impact_severity	biotype	is_exonic	is_coding	is_lof" > exp

gemini query -q "select gene, impact, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variant_impacts" \
                    --header \
                    test1.vep.db \
					| sort -k1,1 -k2,2 > obs
check obs exp
exit
rm obs exp

#########################################################################
# 4. Test variants table for severe_impact columns (VEP)
#########################################################################
echo "    effstring.t04...\c"
echo "anno_id	gene	transcript	impact	impact_severity	biotype	is_exonic	is_coding	is_lof
1	MIR1302-10	ENST00000473358	intron_variant	LOW	lincRNA	0	0	0
1	OR4F5	ENST00000335137	synonymous_variant	LOW	protein_coding	1	1	0
1	OR4F5	ENST00000335137	missense_variant	MED	protein_coding	1	1	0
2	SAMD11	ENST00000342066	frameshift_variant	HIGH	protein_coding	1	1	1
2	SAMD11	ENST00000342066	3_prime_UTR_variant	LOW	protein_coding	1	0	0
2	HES4	ENST00000428771	5_prime_UTR_variant	LOW	protein_coding	1	0	0
2	CCNL2	ENST00000400809	splice_acceptor_variant	HIGH	protein_coding	0	0	1" > exp

gemini query -q "select anno_id, gene, transcript, impact, impact_severity, biotype, \
                    is_exonic, is_coding, is_lof from variants" \
                    --header \
                    test1.vep.db \
                    > obs
check obs exp
rm obs exp

###########################################################################
