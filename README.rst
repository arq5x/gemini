gemini - a lightweight db framework for disease and population genetics.
======================================================================

``gemini`` is not yet ready for safe consumption, but stay tuned...we're working hard on it.

If you are okay living dangerously and potentially being disappointed, you can install and play with ``gemini`` as follows::

    python setup.py install
    sudo ./install-data.py
    


Overview
--------
The vision for gemini is simple: to provide a simple, flexible, and powerful
framework for exploring genetic variation for disease and population genetics.
We aim to utilize the expressive power of databases (particularly SQL) while
attempting to overcome the fundamental challenges associated with using 
databases for very large (e.g. 1,000,000 variants times 1,000 samples 
yields one billion genotypes) datasets.

``Pop`` will support standard population genetics metrics, genetic disease models,
and provide a wealth of useful annotations.  It will allow the user to explore genetic
variation using custom SQL queries, a large set of pre-defined analysis "shortcuts",
and it will eventually support data visualization and interactive plotting.

In short, we intend for it to be a "one stop shop" for large scale human population research.


The workflow will be as follows:

1. Import a VCF file into the ``gemini`` framework. We recommend first annotating your VCF with SnpEff (other tools will be supported soon)::
    
    gemini load -v my.snpEff.vcf my.db
    
2. Explore the variation therein using shortcuts, custom queries, etc.  Here are a few brief examples

- compute the transition / transversion ratio::
  
    gemini get -s tstv my.db
  
- compute the site frequency spectrum::
  
    gemini get -s sfs my.db
  
- extract all transitions with a call rate > 95%::
  
    gemini get -q "select * from variants where sub_type = 'ts' and call_rate >= 0.95" my.db
  
- extract all loss-of-function variants with an alternate allele frequency < 1%::
  
    gemini get -q "select * from variants where is_lof = 1 and aaf >= 0.01" my.db
  
- extract the nucleotide diversity for each variant::
  
    gemini get -q "select chrom, start, end, pi from variants" my.db
  
- combine ``gemini`` with ``bedtools`` to compute nucleotide diversity estimates across 100kb windows::

    gemini get -q "select chrom, start, end, pi from variants" my.db | \
    bedtools map -a hg19.windows.bed -b - -c 4 -o mean


Urgent areas of improvement
---------------------------
1. Full SQL support for the BLOB gts, gt_types, and gt_phases columns.  Currently, some
support is present, but the SQL "parser/interceptor" doesn't handle all cases.  The
goal is to be able to allow "slicing" into the BLOB's to support retrieval of a individual genotypes.  E.g.::

    gemini get -q "select chrom, start, ref, alt, gts.NA12878, gts.NA12879 from variants" my.db
    chr1	100	A	G	A/A	A/G
    chr1	200	C	G	C/G	G/G

2. Support for multiple third-party "functional annotation" tools.  ``gemini`` depends upon tools like ``SnpEff``, 
``VEP``, and ``ANNOVAR`` for predicting the impact of variants on genes.  Currently, we support SnpEff, but our
goal is to support other tools as well.  Currently, this is a bit complex as these tools are changing rapidly, 
and each tool reports functional consequences a bit differently.

3. Add a table that stores a vector of genotypes for each sample.  This will facilitate fast computation of many
popgen metrics.

4. Move variant impacts to a separate table such that the ``variants`` table only contains one row per variant?
