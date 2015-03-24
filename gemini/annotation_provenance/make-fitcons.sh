PRE=hg19_fitcons_fc-i6-0_V1-01
bigWigToBedGraph ${PRE}.bw ${PRE}.bed
bgzip ${PRE}.bed
tabix ${PRE}.bed.gz
