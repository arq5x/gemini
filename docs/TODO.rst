Urgent
======
1. How best to deal with multiple chrom namin schemes?
2. How best to track and handle multiple genome versions.
3. effect ranking system so that group by can be for highest impact
Ideas
=====
4. Popgen:
	a. coefficient of inbreeding
	b. 
5. Visualization
	a. make plots of a given region/genes with snps on top, as well as the allele frequencies, etc.
	b. use genometools libraries.
	c. manhattan plots.
6. Use Pytables for really big files, as the genotype table can be stored as fixed-width records.

7. Interface
	pop load
	pop get
	    # use a raw SQL query
		pop get -q "SELECT..."
		
		# use short cuts that automatically generate the appropriate queries behind the scenes.
		pop get -s variants
		pop get -s tstv
		pop get -s sfs --binsize 0.01
		pop get -s samples
		pop get -s lof
	pop stats
		fst  - model after vcftools?
	
