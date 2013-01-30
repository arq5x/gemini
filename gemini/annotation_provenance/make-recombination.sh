wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz

tar -zxvf genetic_map_HapMapII_GRCh37.tar.gz --exclude=README.txt

# The headers are as follows: Chrom	Position(bp)	Rate(cM/Mb)	Map(cM)
# combine the files and get rid of the headers
# the awk statement is required to create BED format,
# as the original file just has a single position.
# in other words, it takes this:
# chr1	55550	2.981822	0.000000
# chr1	82571	2.082414	0.080572
# chr1	88169	2.081358	0.092229
# chr1	254996	3.354927	0.439456
# chr1	564598	2.887498	1.478148
#
# and creates this:
# chr1	0	55550	2.981822	0.000000
# chr1	55550	82571	2.082414	0.080572
# chr1	82571	88169	2.081358	0.092229
# chr1	88169	254996	3.354927	0.439456
# chr1	254996	564598	2.887498	1.478148
cat genetic_map_GRCh37_* | \
   grep "^chr"  | \
   awk -v start=0 '{OFS="\t"; print $1,start,$2,$3,$4; start=$2}' | \
   sort -k1,1 -k2,2n | \
   bgzip > genetic_map_HapMapII_GRCh37.gz

tabix -p bed genetic_map_HapMapII_GRCh37.gz

# move the annotations to the server and cleanup
scp genetic_map_HapMapII_GRCh37.gz* fir.itc.virginia.edu:~/public_html/files/gemini/annotations/
rm genetic_map_HapMapII_GRCh37*