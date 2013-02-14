###########################################################################################
#1. Test ENCODE tfbs clusters
###########################################################################################
echo "    encode.t1...\c"
echo "None
None
None
E2F6_1,ZNF263_1
None
E2F6_(H-50)_1,Pol2_14,Pol2-4H8_1,HA-E2F1_1,HEY1_2,TBP_3,E2F6_1,TCF4_1,HMGN3_1,TAF1_1,Egr-1_1
GABP_1" > exp
gemini query -q "select encode_tfbs from variants" test1.snpeff.db > obs
check obs exp
rm obs exp

###########################################################################################
#2. Test ENCODE DnaseI clusters
###########################################################################################
echo "    encode.t2...\c"
echo "None	None
None	None
None	None
5	Adult_Th0;Fibrobl;Gm19238;HMEC;Osteobl
6	A549;HAEpiC;HCPEpiC;HRE;HVMF;hESCT0
79	8988t;A549;AG04450;AG09309;AG09319;Adult_Th0;AoAF;AosmcSerumfree;BJ;CD14;Cll;Fibrobl;Fibrop;GM12864;GM12878;Gm18507;Gm19238;Gm19239;Gm19240;H9es;HAEpiC;HAc;HAh;HBMEC;HCF;HCFaa;HCM;HCPEpiC;HCT116;HConF;HEEpiC;HESC;HFF_MyC;HGF;HL60;HMF;HMVEC_dAd;HMVEC_dBlAd;HMVEC_dLyNeo;HMVEC_dNeo;HNPCEpiC;HPAEC;HPAF;HPF;HPdLF;HRE;HRGEC;HSMM;HSMM_D;HUVEC;HVMF;Hepatocytes;Hpde6e6e7;Huh75;Ips;Ishikawa_E;LNCap;Medullo;Melano;Myometr;NHA;NHDF_Ad;NHDF_Neo;NHLF;NT2_D1;Osteobl;Phte;PrEC;Progfib;RPTEC;SAEC;SK_N_MC;SkMC;Stellate;T47d;WI_38_TAM;hESCT0;hTH1;hTH2
7	CMK;Gm19238;HL60;HUVEC;Hela;LNCap;hTH1" > exp
gemini query -q "select encode_dnaseI_cell_count, encode_dnaseI_cell_list from variants" test1.snpeff.db > obs
check obs exp
rm obs exp

############################################################################################
#3. Test encode_consensus
############################################################################################
echo "    encode.t3...\c"
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

######################################################################################
##3. Test encode_segway
######################################################################################
#echo "    encode.t3...\c"
#echo "Low6	Quies	Quies	Repr2	Quies	PromF
#Quies	Quies	Repr2	Repr2	Quies	Low3
#Quies	Quies	Repr1	Quies	Low1	EnhF
#Repr2	PromP	EnhPr	PromP	Repr3	Enh
#Low4	DnaseD	Gen3'	ElonW	Gen3'	Gen3'2
#Enh	Tss	Regul	Regul	Regul	Regul
#TssF	TssF	PromF	PromF	TssF	TssF" > exp
#gemini query -q "select encode_segway_gm12878,encode_segway_h1hesc,encode_segway_helas3, \
#                 encode_segway_hepg2,encode_segway_huvec,encode_segway_k562 \
#                 from variants" test1.snpeff.db > obs
#check obs exp
#rm obs exp                 
#
######################################################################################
##4. Test encode_chromhmm
######################################################################################
#echo "    encode.t4...\c"
#echo "Quies	Quies	Quies	Quies	Quies	Quies
#Quies	Quies	Quies	Quies	Quies	PromF
#Quies	Quies	Quies	Quies	Quies	PromF
#ReprD	ReprD	EnhW	PromP	DnaseD	EnhW
#Gen3'	DnaseD	Pol2	Pol2	DnaseD	Gen5'
#PromP	PromP	Tss	Tss	Tss	Tss
#Tss	Tss	Tss	Tss	Tss	Tss" > exp
#gemini query -q "select encode_chromhmm_gm12878,encode_chromhmm_h1hesc, \
#                encode_chromhmm_helas3,encode_chromhmm_hepg2,encode_chromhmm_huvec, \
#                encode_chromhmm_k562 from variants" test1.snpeff.db > obs
#check obs exp
#rm obs exp