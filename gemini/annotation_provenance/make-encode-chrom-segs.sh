# 1. Get the consensus ENCODE segmentations from EBI.
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/gm12878.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/h1hesc.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/helas3.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/hepg2.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/huvec.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/k562.combined.bb

# 2. Make BEDGRAPHs of the ENCODE segmentation BigBeds
for bigbed in `ls *.combined.bb`
do
   bigBedToBed $bigbed stdout | cut -f 1-4 | bgzip > $bigbed.bedg.gz
done

# 3. Use bedtools to make the union of the ENCODE segmentations.

bedtools unionbedg -i gm12878.combined.bb.bedg.gz \
                      h1hesc.combined.bb.bedg.gz \
                      helas3.combined.bb.bedg.gz \
                      hepg2.combined.bb.bedg.gz \
                      huvec.combined.bb.bedg.gz \
                      k562.combined.bb.bedg.gz \
                    -header \
                    -names gm12878 \
                           h1hesc \
                           helas3 \
                           hepg2 \
                           huvec \
                           k562 \
                    -filler unknown \
          | bgzip \
          > encode.6celltypes.chromseg.bedg.gz

# 4. take a peek (gzcat is for OSX, use zcat elsewhere)
#
#  Glossary for ENCODE chromatin segment predictions.  Taken verbatim from Table 3
# of doi:10.1038/nature11247
#    CTCF:    CTCF-enriched element
#    E:       Predicted enhancer
#    PF:      Predicted promoter flanking region
#    R:       Predicted repressed or low-activity region
#    TSS:     Predicted promoter region including TSS
#    T:       Predicted transcribed region
#    WE:      Predicted weak enhancer or open chromatin cis-regulatory element
#    unknown: added by us.  hopefully self-explanatory
# 
# (gzcat encode.6celltypes.chromseg.bedg.gz  | head -1; gzcat  
# encode.6celltypes.chromseg.bedg.gz | \
#                                                             awk ' > 100000 && NR < 100020')
# chrom	start	end	gm12878	h1hesc	helas3	hepg2	huvec	k562
# chr1	21710800	21711000	WE	R	R	R	R	R
# chr1	21711000	21711200	WE	unknown	R	R	R	R
# chr1	21711200	21711298	WE	unknown	R	R	R	unknown
# chr1	21711298	21711400	WE	unknown	R	R	R	WE
# chr1	21711400	21711468	WE	R	R	R	R	WE
# chr1	21711468	21711779	WE	R	R	R	R	unknown
# chr1	21711779	21711785	unknown	R	R	R	R	unknown
# chr1	21711785	21711900	unknown	R	R	R	R	WE
# chr1	21711900	21712200	unknown	R	R	R	R	unknown
# chr1	21712200	21713000	R	R	R	R	R	R
# chr1	21713000	21713087	R	T	R	R	R	unknown
# chr1	21713087	21713195	R	T	T	R	R	unknown
# chr1	21713195	21713400	R	T	R	R	R	unknown
# chr1	21713400	21714226	R	R	R	R	R	R
# chr1	21714226	21716200	R	R	R	R	T	R
# chr1	21716200	21716205	R	R	R	R	T	unknown
# chr1	21716205	21716288	R	R	R	R	R	unknown
# chr1	21716288	21716425	R	R	R	R	R	E
# chr1	21716425	21716523	T	R	R	R	R	E


# 5. tabix the 6-way consensus map for use within gemini.

tabix -p bed encode.6celltypes.chromseg.bedg.gz