# 1. Get the ENCODE segmentations from EBI.
# consensus
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/gm12878.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/h1hesc.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/helas3.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/hepg2.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/huvec.combined.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/k562.combined.bb

# SegWay
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/gm12878.segway.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/h1hesc.segway.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/helas3.segway.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/hepg2.segway.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/huvec.segway.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/k562.segway.bb

# ChromHMM
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/gm12878.ChromHMM.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/h1hesc.ChromHMM.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/helas3.ChromHMM.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/hepg2.ChromHMM.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/huvec.ChromHMM.bb
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/k562.ChromHMM.bb


# 2. Make BEDGRAPHs of the ENCODE segmentation BigBeds
for bigbed in `ls *.bb`
do
   bigBedToBed $bigbed stdout | cut -f 1-4 | bgzip > $bigbed.bedg.gz
done

# 3. Use bedtools to make the union of each ENCODE segmentation set.
#    That is, make a one file for the consensus segmentations including all 6 cell
#    lines, another for segway, and another for ChromHMM
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
          > encode.6celltypes.consensus.bedg.gz

bedtools unionbedg -i gm12878.segway.bb.bedg.gz \
                    h1hesc.segway.bb.bedg.gz \
                    helas3.segway.bb.bedg.gz \
                    hepg2.segway.bb.bedg.gz \
                    huvec.segway.bb.bedg.gz \
                    k562.segway.bb.bedg.gz \
                  -header \
                  -names gm12878 \
                         h1hesc \
                         helas3 \
                         hepg2 \
                         huvec \
                         k562 \
                  -filler unknown \
        | bgzip \
        > encode.6celltypes.segway.bedg.gz

bedtools unionbedg -i gm12878.ChromHMM.bb.bedg.gz \
                h1hesc.ChromHMM.bb.bedg.gz \
                helas3.ChromHMM.bb.bedg.gz \
                hepg2.ChromHMM.bb.bedg.gz \
                huvec.ChromHMM.bb.bedg.gz \
                k562.ChromHMM.bb.bedg.gz \
              -header \
              -names gm12878 \
                     h1hesc \
                     helas3 \
                     hepg2 \
                     huvec \
                     k562 \
              -filler unknown \
    | bgzip \
    > encode.6celltypes.ChromHMM.bedg.gz

# 5. tabix the 6-way segmentation maps for use within gemini.

tabix -p bed encode.6celltypes.consensus.bedg.gz
tabix -p bed encode.6celltypes.segway.bedg.gz
tabix -p bed encode.6celltypes.chromhmm.bedg.gz
