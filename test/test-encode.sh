###########################################################################################
#1. Test ENCODE tfbs
###########################################################################################
echo "    encode.t1...\c"
echo "None
None
None
E2F6_210_1,ZNF263_334_1
None
E2F6_(H-50)_1000_1,Pol2_1000_14,Pol2-4H8_845_1,HA-E2F1_686_1,HEY1_349_2,TBP_152_3,E2F6_1000_1,TCF4_213_1,HMGN3_434_1,TAF1_247_1,Egr-1_914_1
GABP_26_1" > exp
gemini query -q "select encode_tfbs from variants" test1.snpeff.db > obs
check obs exp
rm obs exp

############################################################################################
#2. Test encode_consensus
############################################################################################
echo "    encode.t2...\c"
echo "R	R	R	R	R	unknown
R	R	R	R	R	unknown
R	R	R	R	R	unknown
R	unknown	WE	E	unknown	WE
T	unknown	T	T	unknown	T
TSS	TSS	TSS	TSS	TSS	TSS
TSS	TSS	unknown	unknown	TSS	TSS" > exp
gemini query -q "select encode_consensus_gm12878,encode_consensus_h1hesc, \
                 encode_consensus_helas3,encode_consensus_hepg2,encode_consensus_huvec, \
                 encode_consensus_k562 from variants" test1.snpeff.db > obs
check obs exp
rm obs exp

#####################################################################################
#3. Test encode_segway
#####################################################################################
echo "    encode.t3...\c"
echo "Low6	Quies	Quies	Repr2	Quies	PromF
Quies	Quies	Repr2	Repr2	Quies	Low3
Quies	Quies	Repr1	Quies	Low1	EnhF
Repr2	PromP	EnhPr	PromP	Repr3	Enh
Low4	DnaseD	Gen3'	ElonW	Gen3'	Gen3'2
Enh	Tss	Regul	Regul	Regul	Regul
TssF	TssF	PromF	PromF	TssF	TssF" > exp
gemini query -q "select encode_segway_gm12878,encode_segway_h1hesc,encode_segway_helas3, \
                 encode_segway_hepg2,encode_segway_huvec,encode_segway_k562 \
                 from variants" test1.snpeff.db > obs
check obs exp
rm obs exp                 

#####################################################################################
#4. Test encode_chromhmm
#####################################################################################
echo "    encode.t4...\c"
echo "Quies	Quies	Quies	Quies	Quies	Quies
Quies	Quies	Quies	Quies	Quies	PromF
Quies	Quies	Quies	Quies	Quies	PromF
ReprD	ReprD	EnhW	PromP	DnaseD	EnhW
Gen3'	DnaseD	Pol2	Pol2	DnaseD	Gen5'
PromP	PromP	Tss	Tss	Tss	Tss
Tss	Tss	Tss	Tss	Tss	Tss" > exp
gemini query -q "select encode_chromhmm_gm12878,encode_chromhmm_h1hesc, \
                encode_chromhmm_helas3,encode_chromhmm_hepg2,encode_chromhmm_huvec, \
                encode_chromhmm_k562 from variants" test1.snpeff.db > obs
check obs exp
rm obs exp