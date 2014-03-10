wget -c http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.4.zip

unzip dbNSFPv2.4.zip

(head -1 dbNSFP2.4_variant.chr1; grep -hv ^# dbNSFP2.4_variant.chr*) \
    | cut -f 1-6,19,20,21,24,26,27,29,30,32,33,35,36,38,39,41,42,44,45,47,48,50 \
    | bgzip \
    > dbNSFP2.4.gemini.txt.gz

tabix -s 1 -b 2 -e 2 dbNSFP2.4.gemini.txt.gz