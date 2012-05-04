gemini - a framework for mining genome variation.
=================================================

``gemini`` is not yet ready for safe consumption, but stay tuned...we're working hard on it.

If you are okay living dangerously and potentially being disappointed, you can install and play with ``gemini`` as follows::

    python setup.py install
    sudo ./install-data.py
    


Overview
--------
The intent of ``gemini`` is to provide a simple, flexible, and powerful
framework for exploring genetic variation for disease and population genetics.
We aim to leverage the expressive power of SQL while attempting to overcome the fundamental challenges associated with using 
databases for very large (e.g. 1,000,000 variants times 1,000 samples 
yields one billion genotypes) datasets.

By augmenting VCF files with many informative annotations, and converting the information
into a ``sqlite`` database framework, ``gemini`` provides a flexible database-driven API for:

- *Ad hoc* queries.  Return all loss-of-function INDELs::

	$ gemini get -q "select * from variants where type = 'indel' and is_lof = 1" my.db

- Pre-defined analysis short-cuts. Compute the ratio of transitions to transversions::

	$ gemini get -s tstv my.db
	transitions	transversions	ts/tv
	1,302,778	511,578		2.547

- A framework for exploring genetic variation in the context of built-in genome annotations. Return the coordinates, alleles, and clinical significance of all OMIM variants with an alternate allele frequency l.t.e 1%::
	
	$ gemini get -q "select chrom, start, ref, alt, clin_sigs \
                     from variants where in_omim = 1 and aaf < 0.1" my.db

- A platform for genomic discovery

- A simple to use framework for developers to create new tools.

Current annotations that are derived for each variant are listed below.  We emphasize that this is a very preliminary list.
Using the framework we have developed, it is very easy for us to add new annotations. Therefore, we intend to add many other annotations
including: ENCODE regulatory and histone modification tracks, conservation, GWAS analyses, etc.:

- cyto_band: the chromosome band based on Giemsa staining
- dbSNP status: is the variant in dbSNP? what are the rsIds?
- OMIM status
- Clinical significance
- RepeatMasker annotations
- Overlap with CpG Islands
- Overlap with segmental duplications



Basic workflow
---------------

1. Import a VCF file into the ``gemini`` framework. We recommend first annotating your VCF with ``SnpEff`` or ``VEP`` (other tools may be supported soon).  
   In the process of loading the VCF into the database framework, many other annotations are calculated for each variant and stored for 
   subsequent querying/analysis::
    
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


Active areas of improvement
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

4. Add ESP allele frequency variation as an annotation!

	- ESP5400.snps.vcf.tar.gz from the ESP variant server website.

6. Add GERP score.

7. Add KEGG

8. Add COSMIC
