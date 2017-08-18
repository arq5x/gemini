wget -O - https://data.broadinstitute.org/gnomAD/release-170228/exomes/vcf/gnomad.exomes.r2.0.1.sites.vcf.gz | zcat - | perl -pe 's/;CSQ.+//g' | less | bgzip -c > gnomad.exomes.r2.0.1.sites.no-VEP.vcf.gz

zcat gnomad.exomes.r2.0.1.sites.no-VEP.vcf.gz \
    | perl -pe 's/(;[^_\s;]+_HIST(_AL[LT])?=[^;]+)//g' \
    | vt decompose -s - \
    | vt normalize -r ~/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa - \
    | bgzip -c > gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz

tabix -f gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz
