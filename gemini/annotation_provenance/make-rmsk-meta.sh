curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz | \
     gzcat | \
     awk '{OFS="\t"; print $6,$7,$8,$12"_"$13"_"$11}' \
     > rmsk.bed

curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz | \
     gzcat | \
     cut -f 2-5 \
     > simrep.bed

curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/microsat.txt.gz | \
     gzcat | \
     cut -f 2-5 \
         > microsat.bed

cat rmsk.bed simrep.bed microsat.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -nms -i - | \
    bgzip \
    > hg19.rmsk.bed.gz

tabix -p bed hg19.rmsk.bed.gz