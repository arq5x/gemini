wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz \
    | zcat - \
    | perl -pe 's/;CSQ.+//g' \
    | perl -pe 's/(;[^_\s;]+_HIST(_AL[LT])?=[^;]+)//g' \
    | ~/src/vt/vt decompose -s - \
    | ~/src/vt/vt normalize -r /data/human/g1k_v37_decoy.fa - \
    | bgzip -c > gnomad.exomes.r2.0.2.sites.no-VEP.nohist.tidy.vcf.gz

tabix -f gnomad.exomes.r2.0.2.sites.no-VEP.nohist.tidy.vcf.gz
