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

# setup the testing databases from the testing VCF files
sh data-setup.sh

# Test basic functionality
sh test-columns.sh

# Test loading functionality
sh test-load.sh

# Test genotype BLOB functionality
sh test-genotypes.sh

# Test ClinVar attributes
sh test-clinvar.sh

# Test population_gen metrics
sh test-pop.sh

# Test mappability
sh test-map.sh

# Test genome annotations
sh test-genome.sh

# Test encode annotations
sh test-encode.sh

# Test EFF string derived elements in INFO column
sh test-effstring.sh

# Test loading functionality
sh test-annotate-tool.sh

# Test pathway tool
sh test-pathtool.sh

# Test interaction tool
sh test-interactions.sh

# cleanup
rm ./*.db

