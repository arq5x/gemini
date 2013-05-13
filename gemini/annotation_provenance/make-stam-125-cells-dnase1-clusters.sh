# data setup
wget ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/combined_peaks/multi-tissue.master.ntypes.simple.hg19.bed
wget http://www.uwencode.org/proj/Nature_Thurman_et_al/data/multi-tissue.master.ntypes.simple.names.hg19.txt.gz
gunzip multi-tissue.master.ntypes.simple.names.hg19.txt.gz


# peak at the sites (col 4 = number of cell types, col 5 = max signal in any one cell type)
# head -5 multi-tissue.master.ntypes.simple.hg19.bed
# chr1	10120	10270	MCV-37	39.232800
# chr1	10440	10590	MCV-4	36.369500
# chr1	16140	16290	MCV-5	18.768900
# chr1	20060	20210	MCV-1	9.082430
# chr1	56200	56350	MCV-1	18.726500

# peak at the list of cell types that were observed to have DnaseI hypersensitivty.  One line per interval in file above
head -5 
# A549;AG04450;CD20;CD34;CMK;GM12864;GM12865;HAEpiC;HAc;HBMEC;HCF;HCFaa;HESC;HL60;HMVEC_LBl;HMVEC_LLy;HMVEC_dBlAd;HMVEC_dBlNeo;HMVEC_dLyAd;HMVEC_dLyNeo;HMVEC_dNeo;HPAEC;HRPEpiC;HSMM_D;HUVEC;HVMF;Hela;HepG2;Jurkat;LNCap;NB4;PANC1;SK_N_MC;WERI_Rb1;hESCT0;hTH1;hTH2
# BE2_C;HESC;Jurkat;NB4
# CD34;CMK;Fibrobl;HFF;SK_N_MC
# SK_N_MC
# HEEpiC

# add the names of the cell types to the count
paste multi-tissue.master.ntypes.simple.* | sed -e 's/MCV\-//' | sort -k1,1 -k2,2n > stam.125cells.dnaseI.hg19.bed
bgzip stam.125cells.dnaseI.hg19.bed
tabix -p bed stam.125cells.dnaseI.hg19.bed.gz