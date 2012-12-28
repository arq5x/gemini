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

# Test genotype BLOB functionality
sh test-genotypes.sh


