echo "vcf output.t1"
gemini query -q "select chrom, start, end  from variants limit 1" test.query.db --format vcf > vcf

echo "1" > exp
grep -c "AC=6;" vcf > obs
check obs exp

echo "vcf output.t1"
grep -c "DP=7;" vcf > obs
check obs exp

rm vcf obs exp
