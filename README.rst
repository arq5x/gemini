pop - a unique database framework for disease and population genetics.
======================================================================

``pop`` is not yet ready for safe consumption, but stay tuned...we're working hard on it.
The vision is outlined below.

Overview
--------
The vision for pop is simple: to provide a simple, flexible, and powerful
framework for exploring genetic variation for disease and population genetics.
We aim to utilize the expressive power of databases (particularly SQL) while
attempting to overcome the fundamental challenges associated with using 
databases for very large (e.g. 1,000,000 variants times 1,000 samples 
yields one billion genotypes) datasets.

``Pop`` will support standard population genetics metrics, genetic disease models,
and provide a wealth of useful annotations.  It will allow the user to explore genetic
variation using custom SQL queries, a large set of pre-defined analysis "shortcuts",
and it will eventually support data visualization and interactive plotting.


The workflow will be as follows:

1. Import a VCF file into the ``pop`` framework.

    pop load -v my.vcf my.db
    
2. Explore the variation therein using shortcuts, custom queries, etc.  Here are a few brief examples:

    # compute the transition / transversion ratio
        pop get -s tstv my.db
    
    # compute the site frequency spectrum
        pop get -s sfs my.db
        
    # extract all transitions with a call rate > 95%
        pop get -q "select * from variants where sub_type = 'ts' and call_rate >= 0.95"
        
    # extract the nucleotide diversity for each variant
        pop get -q "select chrom, start, end, pi from variants"
        
    # combine ``pop`` with ``bedtools`` to compute nucleotide diversity estimates across 100kb windows
        pop get -q "select chrom, start, end, pi from variants" | \
        bedtools map -a hg19.windows.bed -b - -c 4 -o mean