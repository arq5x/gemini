#/usr/bin/sh

# grab the ENCODE TFBS clusters from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz

# add a column describing HOW MANY cells had ChIPseq peaks for the given TF
# column 15 in the BED file is a list of the experiment "scores" for each of the
# cells tested.  Here, we are just pre-computing the count of cells that had non-zero scores.
# Further information on this file can be found at the UCSC Table Browser for table: wgEncodeRegTfbsClusteredV2
gunzip wgEncodeBroadHmmGm12878HMM.bed.gz
cat wgEncodeRegTfbsClusteredV2.bed | awk '
{ 
  OFS="\t"
  split($15, cell_scores, ",")
  cell_count = 0
  for (i=0; i < length(cell_scores); i++) {
    if (cell_scores[i] > 0) {
      cell_count += 1
    }
  }
  print $1,$2,$3,$4"_"$5"_"cell_count,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15
}' > wgEncodeRegTfbsClusteredV2.cell_count.bed

# tabix
bgzip wgEncodeRegTfbsClusteredV2.cell_count.bed
tabix -p bed wgEncodeRegTfbsClusteredV2.cell_count.bed.gz