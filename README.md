gemini - a framework for mining genome variation.
=================================================

``gemini`` is not yet ready for safe consumption, but stay tuned...we're working hard on it.

If you are okay living dangerously and potentially being disappointed, you can install ``gemini`` as follows:

    python setup.py install

One of the more appealing features in ``gemini`` is that it automatically annotates variants in a VCF file with several
genome annotations (see below for more details).  However, you must first install these data files on your system. 
It's easy enough --- you just need to run the following script and tell it in which what full path you'd like to install
the necessary data files. The recommended path is in ``/usr/local/share``:

    gemini/install-data.py /usr/local/share/gemini


Overview
========
The intent of ``gemini`` is to provide a simple, flexible, and powerful
framework for exploring genetic variation for disease and population genetics.
We aim to leverage the expressive power of SQL while attempting to overcome the fundamental challenges associated with using 
databases for very large (e.g. 1,000,000 variants times 1,000 samples 
yields one billion genotypes) datasets.

Have a look at this [poster](http://dl.dropbox.com/u/515640/posters_and_slides/Quinlan-Gemini-Poster.pdf) to get a high-level sense of what ``gemini`` is trying to accomplish.

The first thing we have to do is load an existing VCF file into the gemini framework.  We expect you to have
annotated the functional consequence of each variant in your VCF using either VEP or snpEff (>= v3.0).  Logically,
the loading step is done with the ``gemini load`` command.  Below are two examples based on a VCF file that
we ingeniously name my.vcf.  The first example assumes that the VCF has been pre-annotated with VEP and the second
assumes snpEff.

	# VEP-annotated VCF
	$ gemini load -v my.vcf -t VEP my.db

	# snpEff-annotated VCF
	$ gemini load -v my.vcf -t snpEff my.db

As each variant is loaded into the ``gemini`` database framework, it is being compared against several
annotation files that come installed with the software.  We have developed an annotation framework
that leverages Tabix, BEDTools, and pybedtools to make things easy and fairly performant. The idea is that, 
by augmenting VCF files with many informative annotations, and converting the information
into a ``sqlite`` database framework, ``gemini`` provides a flexible database-driven API for data exploration,
visualization, population genomics and medical genomics.  **We feel that this ability to integrate variation
with the growing wealth of gemome annotations is the most compelling aspect of gemini ``gemini``**.  Combining this
with the ability to explore data with SQL using a database design that can scale to 1000s of individuals (genotypes too!)
makes for a tasty data exploration cupcake.  Here are some examples of the things you can do.

*Ad hoc* queries.  Return all loss-of-function INDELs.

    $ gemini query -q "select * from variants where type = 'indel' and is_lof = 1" my.db

Pre-defined analysis short-cuts. Compute the ratio of transitions to transversions::

    $ gemini stats --tstv my.db
    transitions	transversions	ts/tv
    1,302,778	511,578		2.547

A framework for exploring genetic variation in the context of built-in genome annotations. Return the coordinates, alleles, and clinical significance of all OMIM variants with an alternate allele frequency l.t.e 1%::
	
    $ gemini query -q "select chrom, start, ref, alt, clin_sigs \
                     from variants where in_omim = 1 and aaf < 0.1" my.db

A platform for genomic discovery
A simple to use framework for developers to create new tools.


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
==============

**Import a VCF file into the ``gemini`` framework.**
We recommend first annotating your VCF with ``SnpEff`` or ``VEP`` (other tools may be supported soon).  In the process of loading the VCF into the database framework, many other annotations are calculated for each variant and stored for subsequent querying/analysis. **Note**: If using snpEff, ee currently require VCFs to be annotated with version 3.0 or later.
    
    gemini load -v my.snpEff.vcf -t snpEff my.db
    
**Explore variation using shortcuts. Here are a few brief examples**

Compute the transition / transversion ratio:
  
    gemini stats --tstv my.db
  
Compute the site frequency spectrum:
  
    gemini stats --sfs my.db

Compute the pairwise genetic distance for use with PCA:

    gemini stats --mds my.db


**Explore variation  using custom queries. Here are a few brief examples:**

Extract all transitions with a call rate > 95%::
  
    gemini query -q "select * from variants where sub_type = 'ts' and call_rate >= 0.95" my.db
  
Extract all loss-of-function variants with an alternate allele frequency < 1%:
  
    gemini query -q "select * from variants where is_lof = 1 and aaf >= 0.01" my.db
  
Extract the nucleotide diversity for each variant:
  
    gemini query -q "select chrom, start, end, pi from variants" my.db

Combine ``gemini`` with ``bedtools`` to compute nucleotide diversity estimates across 100kb windows:

    gemini query -q "select chrom, start, end, pi from variants order by chrom, start, end" my.db | \
    bedtools map -a hg19.windows.bed -b - -c 4 -o mean

Adding your own, custom annotations to the ``gemini`` framework
===========================================================
It is inevitable that researchers will want to enhance the ``gemini``
framework with their own, custom annotations. ``gemini`` provides a
sub-command called ``annotate`` for exactly this purpose.  As long as 
you provide a ``tabix``'ed annotation file in either BED or VCF format,
the ``annotate`` tool will, for each variant in the ``variants`` table,
screen for overlaps in your annotation file and update a new column in the
``variants`` table that you may specify on the command line.  This is best
illustrated by example. 

Let's assume you have already created a ``gemini`` database of a VCF file
using the ``load`` module.

	gemini load -v my.vcf -t VEP my.db

Add a new column called "my_col" that tracks whether a given variant overlapped (1)
or did not overlap (0) intervals in your annotation file.

	# Add my_col to the database as a _boolean_
	gemini annotate -f my_annos.bed.gz -c my_col -t boolean my.db

	# Now query the results
	gemini query -q "select chrom, start, end, variant_id, my_col from variants" my.db | head -3
	chr22	16504479	16504480	1	1
	chr22	16504488	16504489	2	0
	chr22	16504490	16504491	3	1

Add a new column called "my_col" that counts the number of overlaps a given variant 
has with intervals in your annotation file.

	# Add my_col to the database as a _count_
	gemini annotate -f my_annos.bed.gz -c my_col -t count my.db

	# Now query the results
	gemini query -q "select chrom, start, end, variant_id, my_col from variants" my.db | head -3
	chr22	16504479	16504480	1	2
	chr22	16504488	16504489	2	0
	chr22	16504490	16504491	3	1

Add a new column called "my_col" that creates a list of a specific column from the
annotation file for each given variant.

	# Add my_col to the database as a _list_, extracting the 4th
	# column from the annotation file to create the list
	gemini annotate -f my_annos.bed.gz -c my_col -t list -e 4 my.db

	# Now query the results
	gemini query -q "select chrom, start, end, variant_id, my_col from variants" my.db | head -3
	chr22	16504479	16504480	1	rs123,rs456
	chr22	16504488	16504489	2	None
	chr22	16504490	16504491	3	rs789


Extracting KEGG pathways for variants/genes/individuals
=======================================================

``gemini`` catalogs KEGG pathway information and using the 
``pathways`` tool, one can extract pathway information for each
sample that has variants affecting a given gene.  The only requirement
is that we know what version of Ensembl genes were using by snpEff or VEP.
Currently, we expect version 66, 67, or 68.  

    gemini pathways -v 66 chr22.low.exome.snpeff.100samples.vcf.db

We can also focus solely on loss-of-function mutations with the ``--lof``
argument:

    gemini pathways --lof -v 66 chr22.low.exome.snpeff.100samples.vcf.db | head
    chrom	start	end	ref	alt	highest_impact	sample	genotype	gene	transcript	pathway
    chr22	18912676	18912677	C	T	stop_gain	HG00312	C|T	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	HG01069	C|T	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	NA12275	C|T	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	NA18535	C|T	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	NA19324	T|C	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	NA19327	T|C	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	NA19655	C|T	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	18912676	18912677	C	T	stop_gain	NA20341	C|T	PRODH	ENST00000357068	hsa01100:Metabolic_pathways,hsa00330:Arginine_and_proline_metabolism
    chr22	19258043	19258044	A	<DEL>	exon_deleted	NA18999	A|<DEL>	CLTCL1	ENST00000505027	hsa04721:Synaptic_vesicle_cycle,hsa04961:Endocrine_and_other_factor_regulated_calcium_reabsorption,hsa05016:Huntington's_disease,hsa05100:Bacterial_invasion_of_epithelial_cells,hsa04142:Lysosome,hsa04144:Endocytosis


Extracting variants from specific regions or genes
==================================================

``gemini`` allows one to extract variants that fall within specific genomic coordinates as follows

	gemini region --reg chr1:100-200 my.db

Or, one can extract variants based on a specific gene name.

	gemini region --gene PTPN22 my.db


Conducting analyses on genome "windows".
========================================

``gemini`` includes a convenient tool for computing variation metrics across genomic windows (fixed and sliding).
Here are a few examples to whet your appetite.  If you're still hungry, email us.

Compute the average nucleotide diversity for all variants found in non-overlapping, 50Kb windows.

	gemini windower -w 50000 -s 0 -t nucl_div -o mean my.db

Compute the average nucleotide diversity for all variants found in 50Kb windows that overlap by 10kb.

	gemini windower -w 50000 -s 10000 -t nucl_div -o mean my.db



Working with individual genotypes.
==================================
This will take some time to refine and explain.  I am being lazy.  To do.


Active areas of improvement
===========================

Full SQL support for the BLOB gts, gt_types, and gt_phases columns.  Currently, some
support is present, but the SQL "parser/interceptor" doesn't handle all cases.  The
goal is to be able to allow "slicing" into the BLOB's to support retrieval of a individual genotypes.  E.g.::

    gemini query -q "select chrom, start, ref, alt, gts.NA12878, gts.NA12879 from variants" my.db
    chr1	100	A	G	A/A	A/G
    chr1	200	C	G	C/G	G/G

Support for multiple third-party "functional annotation" tools.  ``gemini`` depends upon tools like ``SnpEff``, 
``VEP``, and ``ANNOVAR`` for predicting the impact of variants on genes.  Currently, we support SnpEff and VEP, but our
goal is to support other tools (such as the VAT from the Gerstein lab) as well.  Currently, this is a bit complex 
as these tools are changing rapidly, and each tool reports functional consequences a bit differently.

Add a table that stores a vector of genotypes for each sample.  This will facilitate fast computation of many
popgen metrics.

Add GERP score.
Add KEGG
Add COSMIC
