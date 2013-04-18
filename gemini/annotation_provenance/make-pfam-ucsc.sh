##############################################################################################
# make input: "Pfam A domains with chromosomal co-ordinates from UCSC table browser"
# Dated: 17April 2013
# -> Go to UCSC table browser:
# -> choose; assembly: hg19, genome: Human, region: genome
# -> Select 'Genes and Gene Prediction Tracks' under group and 'Pfam in UCSC Gene' under track
# -> select 'BED- browser extensible data' under output format and click 'get output'
# -> In the resulting page, select 'Exons plus' under 'Create one BED record per:' and 
# -> click get BED
##############################################################################################

cut -f1-4 input > hg19.pfam.ucscgenes
sed 's/\_exon.*$//g' hg19.pfam.ucscgenes > hg19.pfam.ucscgenes.bed
cat hg19.pfam.ucscgenes.bed | sort -k1,1 -k2,2n | uniq | bgzip -c > hg19.pfam.ucscgenes.bed.gz
tabix -p bed hg19.pfam.ucscgenes.bed.gz
