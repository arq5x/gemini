
echo "test.multi-col.t1"

echo "variant_id	chrom	start	ref	alt	gene	biotype	gene	variant_id	biotype	transcript
1	chr1	30547	T	G	FAM138A	protein_coding	FAM138A	1	protein_coding	ENST00000417324
1	chr1	30547	T	G	FAM138A	protein_coding	FAM138A	1	protein_coding	ENST00000461467" > exp

gemini query --header -q "select v.variant_id, v.chrom, v.start, v.ref, v.alt, v.gene, v.biotype, i.gene, i.variant_id, i.biotype, i.transcript from variants\
                          v, variant_impacts i where v.variant_id=i.variant_id limit 2" test.query.db > obs
check obs exp
