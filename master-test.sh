check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

cd test
rm ./*.db

# setup the testing databases from the testing VCF files
bash data-setup.sh

bash test-mendel-error.sh

bash test-genotype-likelihoods.sh

# Test gemini region
bash test-region.sh

# Test amending the database
bash test-amend.sh

# Test query tool
bash test-query.sh

# Test database dumping
bash test-dump.sh

# Test burden tests
bash test-burden.sh

# Test basic functionality
bash test-columns.sh

# Test cadd scores
bash test-cadd.sh

# Test cosmic
bash test-cosmic.sh

# Test loading functionality
bash test-load.sh

# Test genotype BLOB functionality
bash test-genotypes.sh

# Test ClinVar attributes
bash test-clinvar.sh

# Test Exac
bash test-exac.sh

# Test population_gen metrics
bash test-pop.sh

# Test mappability
bash test-map.sh

# Test genome annotations
bash test-genome.sh

# Test encode annotations
bash test-encode.sh

# Test EFF string derived elements in INFO column
bash test-effstring.sh

# Test loading functionality
bash test-annotate-tool.sh

# Test comp_hets tool
bash test-comphet.sh

# Test pathway tool
bash test-pathtool.sh

# Test interaction tool
bash test-interactions.sh

# Test lof sieve tool
bash test-lofsieve.sh

# Test stats tool
bash test-stats.sh

# Test windower
bash test-windower.sh

#
bash test-fitcons.sh

# Test pfam domains
bash test-pfam.sh

# Test GERP scores
bash test-gerp.sh

# Test disease models
bash test-auto-dom.sh
bash test-auto-rec.sh

# Test denovo tool
bash test-de-novo.sh

# Test wildcards
bash test-wildcards.sh

# Test ROH
bash test-roh.sh

# Test somatic variants
bash test-somatic.sh

# Test fusions
bash test-fusions.sh

bash test-multiple-alts.sh

bash test-bcolz.sh

bash test-esp.sh

bash test-dashes.sh

# backwards compat with no PL/GL columns
bash test-no-gls.sh

bash test-multi-col.sh

# cleanup
rm ./*.db
rm ./file.dot

