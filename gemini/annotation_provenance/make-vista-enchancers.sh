# extract positive human enhancers from VISTA
curl -s "http://enhancer.lbl.gov/cgi-bin/imagedb3.pl?show=1;page=1;form=ext_search;search.status=Positives;search.org=Human;page_size=100;search.result=yes;action=search;search.gene=;search.sequence=1" \
  | grep "^>" \
  > vista.human.positive.enhancers.headers.fa

# convert the VISTA FASTA headers to BED
cat vista.human.positive.enhancers.headers.fa \
  | tr -d "|" \
  | perl -ne 's/\s+/\t/g; print "$_\n" ' \
  | sed -e 's/\>Human//' \
  | tr ":" "\t" \
  | tr "-" "\t" \
  | awk -v OFS="\t" '{list="positives:"; for(i=7; i<=NF; ++i) {list=list$i","}; print $1,$2,$3,$4"_"$5,list}' \
  | sort -k1,1 -k2,2n \
  > hg19.vista.enhancers.20131108.bed

# TABIX
bgzip hg19.vista.enhancers.20131108.bed
tabix -p bed hg19.vista.enhancers.20131108.bed.gz