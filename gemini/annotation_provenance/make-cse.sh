# CSE: Context specific error regions
# https://code.google.com/p/discovering-cse/

wget https://discovering-cse.googlecode.com/files/cse-tracks.tar.gz
tar -xzvpf cse-tracks.tar.gz
cd cse-tracks
mv HiSeq-hg-8_4.bed cse-hiseq-8_4-2013-02-20.bed
bgzip cse-hiseq-8_4-2013-02-20.bed
tabix -p bed cse-hiseq-8_4-2013-02-20.bed.gz
