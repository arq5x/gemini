set -e
wget -O - -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
     | zcat \
     | awk '{FS=OFS="\t"; print $6,$7,$8,$12"_"$13}' \
     | python convert-rmsk.py \
     > rmsk.bed

wget -O - -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz \
     | zcat \
     | awk 'BEGIN{FS=OFS="\t"}{ print $2,$3,$4,"simple"}' \
     > simrep.bed

wget -O - -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/microsat.txt.gz \
     | zcat \
     | awk 'BEGIN{FS=OFS="\t"}{ print $2,$3,$4,"microsat"}' \
     > microsat.bed

wget -O - -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz \
    | zcat \
    | awk 'BEGIN{FS=OFS="\t"}{ print $2,$3,$4,"segdup\n"$8,$9,$10,"segdup" }' \
    > segdup.bed

cat rmsk.bed simrep.bed microsat.bed segdup.bed \
    | perl -pe 's/^chr//' \
    | sort -u -k1,1 -k2,2n  \
    | bgzip \
    > repeats.37.bed.gz
tabix -p bed repeats.37.bed.gz

zcat repeats.37.bed.gz | cut -f 4 | perl -pe 's/,/\n/g' | sort | uniq -c
