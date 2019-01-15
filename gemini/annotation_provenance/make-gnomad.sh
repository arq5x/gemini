set -euo pipefail
<<DONE
wget --quiet https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz

bcftools view --threads 3 gnomad.exomes.r2.1.sites.vcf.bgz -O b -o gnomad.exomes.r2.1.sites.bcf
wait

rm gnomad.exomes.r2.1.sites.vcf.bgz
bcftools index --threads 3 gnomad.exoms.r2.1.sites.bcf
DONE

bcftools annotate -x INFO/variant_type,INFO/allele_type,QUAL,ID,INFO/vep,$(bcftools view -h gnomad.exomes.r2.1.sites.bcf | grep -Po "ID=.*?(hist|RankSum|faf|topmed)[^\,]*" | awk '{ print substr($0, 4, length($0)) }' | perl -pe 's/(.+)\n/INFO\/$1,/g' | perl -pe 's/,$//') \
  --threads 5 -O z -o gnomad.exomes.r2.1.tidy.bcf -O b gnomad.exomes.r2.1.sites.bcf
bcftools index --threads 5 -m 10 gnomad.exomes.r2.1.tidy.bcf
