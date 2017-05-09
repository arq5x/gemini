#!/bin/bash
# Download the new, and final Roadmap ChromHMM data, more specifically, the 15 chromatin states model (primary chromHMM),for 127 tissue and cell types
# All the files are compressed in one file 
wget -c http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz
# Uncompress
tar -zxvf all.mnemonics.bedFiles.tgz
#sort all the bed files using fast sort from bedops
 for i in E*.bed;do sort-bed $i >sorted_$i
 #Optionally we can compress all the sorted files
 ls sorted*|xargs -n1 bgzip
#Merge all the files in one bed file with bedtools unionbedg giving them a header line with their Roadmap IDS 
#instead of their long standard names, this will considerably shorten the query
 bedtools unionbedg -i sorted_E* -header -names E{001..127} >All_chromHMM.bed
# This step will create a large file 11G! so make sure you have enough space
# Now the final compression
bgzip All_chromHMM.bed
#Amazingly the final file is only 294MB
tabix -p bed All_chromHMM.bed.gz




 