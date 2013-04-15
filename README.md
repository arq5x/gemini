GEMINI - a framework for mining genome variation.
=================================================

Overview
========
The intent of ``GEMINI`` is to provide a simple, flexible, and powerful
framework for exploring genetic variation for disease and population genetics.
We aim to leverage the expressive power of SQL while attempting to overcome
the fundamental challenges associated with using databases for very large
(e.g. 1,000,000 variants times 1,000 samples yields one billion genotypes)
datasets.

**NOTE:**  GEMINI is currently in a beta state as we move towards publication.
We welcome comments and suggestions.

Documentation
================

The official documentation is here: http://gemini.readthedocs.org/en/latest/

Acknowledgements
================
GEMINI is being developed in the Quinlan lab (quinlanlab.org) at the University
of Virginia and is led by Uma Paila.  Substantial contributions have also been
made by Brad Chapman (@chapmanb), Rory Kirchner (@roryk), and Oliver Hofmann
at the Harvard School of Public Health.

Installation
============
Install ``GEMINI`` using the automated installation script::

    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
    python gemini_install.py /usr/local /usr/local/share/gemini

This installs GEMINI along with required python libraries, third party tools and
data files used for variant annotation. The installation documentation contains
additional details on installed files and tools.

GEMINI Browser
================
Currently, the majority of GEMINI's functionality is available via a command-line
interface.  However, we are developing a browser-based interface for easier exploration
of GEMINI databases created with the `gemini load` command.

Ironically, as of now, one must launch said browser from the command line as
follows (where `my.db` should be replaced with the name of the GEMINI database
you would like to explore).

    $ gemini browser my.db

At this point, the GEMINI browser is running on port 8088 on your local machine.
Open a web browser to http://localhost:8088/.  You should see something like:

![](https://raw.github.com/arq5x/gemini/master/img/browser-query.png)


Caveats and Limitations
=======================
``GEMINI`` is currently in an *alpha* state.  The basic functionality that we
intend to develop is present, but we may tweak the API and/or restructure
the underlying database as needed.






Have a look at this [poster](http://dl.dropbox.com/u/515640/posters_and_slides/Quinlan-Gemini-Poster.pdf) to get a high-level sense of what ``GEMINI`` is trying to accomplish.

The first thing we have to do is load an existing VCF file into the gemini framework.  We expect you to have
annotated the functional consequence of each variant in your VCF using either VEP or snpEff (Note that v3.0+ of snpEff is required to track the amino acid
length of each impacted transcript).  Logically,
the loading step is done with the ``gemini load`` command.  Below are two examples based on a VCF file that
we ingeniously name my.vcf.  The first example assumes that the VCF has been pre-annotated with VEP and the second
assumes snpEff.

	# VEP-annotated VCF
	$ gemini load -v my.vcf -t VEP my.db

	# snpEff-annotated VCF
	$ gemini load -v my.vcf -t snpEff my.db

Now, the loading step is very computationally intensive and thus can be very slow
with just a single core.  However, if you have more CPUs in your arsenal,
you specify more cores.  This provides a roughly linear increase in speed as a
function of the number of cores. **NOTE: using multiple cores requires that
you have both the `bgzip` tool from `tabix` and the `grabix`
(https://github.com/arq5x/grabix) tool installed in your PATH**.  On our local
machine, we are able to load a VCF file derived from the exomes of 60 samples
in about 10 minutes.  With a single core, it takes a few hours.

	gemini load -v my.vcf -t snpEff --cores 20 my.db

As each variant is loaded into the ``GEMINI`` database framework, it is being compared against several
annotation files that come installed with the software.  We have developed an annotation framework
that leverages Tabix, BEDTools, and pybedtools to make things easy and fairly performant. The idea is that,
by augmenting VCF files with many informative annotations, and converting the information
into a ``sqlite`` database framework, ``GEMINI`` provides a flexible database-driven API for data exploration,
visualization, population genomics and medical genomics.  **We feel that this ability to integrate variation
with the growing wealth of gemome annotations is the most compelling aspect of GEMINI ``GEMINI``**.  Combining this
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

**Import a VCF file into the ``GEMINI`` framework.**
We recommend first annotating your VCF with ``SnpEff`` or ``VEP`` (other tools may be supported soon).  In the process of loading the VCF into the database framework, many other annotations are calculated for each variant and stored for subsequent querying/analysis. **Note**: If using snpEff, ee currently require VCFs to be annotated with version 3.0 or later.

    $ gemini load -v my.snpEff.vcf -t snpEff my.db

**Explore variation using shortcuts. Here are a few brief examples**

Compute the transition / transversion ratio:

    $ gemini stats --tstv my.db

Compute the site frequency spectrum:

    $ gemini stats --sfs my.db

Compute the pairwise genetic distance for use with PCA:

    $ gemini stats --mds my.db


**Explore variation  using custom queries. Here are a few brief examples:**

Extract all transitions with a call rate > 95%::

    $ gemini query -q "select * from variants where sub_type = 'ts' and call_rate >= 0.95" my.db

Extract all loss-of-function variants with an alternate allele frequency < 1%:

    $ gemini query -q "select * from variants where is_lof = 1 and aaf >= 0.01" my.db

Extract the nucleotide diversity for each variant:

    $ gemini query -q "select chrom, start, end, pi from variants" my.db

Combine ``GEMINI`` with ``bedtools`` to compute nucleotide diversity estimates across 100kb windows:

    $ gemini query -q "select chrom, start, end, pi from variants order by chrom, start, end" my.db | \
    bedtools map -a hg19.windows.bed -b - -c 4 -o mean


Queries selecting and filtering upon individual genotypes
===========================================================
The above examples illustrate *ad hoc* queries that do not request or filter
upon the genotypes of individual samples.  Since ``GEMINI`` stores the genotype
information for each variant in compressed arrays that are stored as BLOBs
in the database, standard SQL queries cannot directly access individual
genotypes. However, we have enhanced the SQL syntax to support such queries
with C "struct-like" access.  For example, to retrieve the alleles for a given
sample's (in this case, sample 1094PC0009), one would add `gts.1094PC0009`
to the select statement.

Here is an example of selecting the genotype alleles for four
different samples (note the examples below use the test.snpEff.vcf.db
file that is created in the ./test directory when you run the
`sh master-test.sh` command as described above)::

    $ gemini query -q "select chrom, start, end, ref, alt, gene, \
	                 gts.1094PC0005, \
					 gts.1094PC0009, \
                     gts.1094PC0012, \
					 gts.1094PC0013 \
                 from variants" test.snpEff.vcf.db
    chr1	30547	30548	T	G	FAM138A	./.	./.	./.	./.
    chr1	30859	30860	G	C	FAM138A	G/G	G/G	G/G	G/G
    chr1	30866	30869	CCT	C	FAM138A	CCT/CCT	CCT/CCT	CCT/C	CCT/CCT
    chr1	30894	30895	T	C	FAM138A	T/C	T/C	T/T	T/T
    chr1	30922	30923	G	T	FAM138A	./.	./.	./.	./.
    chr1	69269	69270	A	G	OR4F5	./.	./.	G/G	G/G
    chr1	69427	69428	T	G	OR4F5	T/T	T/T	T/T	T/T
    chr1	69510	69511	A	G	OR4F5	./.	./.	A/G	A/G
    chr1	69760	69761	A	T	OR4F5	A/A	A/T	A/A	A/A
    chr1	69870	69871	G	A	OR4F5	./.	G/G	G/G	G/G

You can also add a header so that you can keep track of who's who::

    $ gemini query -q "select chrom, start, end, ref, alt, gene,
                     gts.1094PC0005, \
    				 gts.1094PC0009, \
                     gts.1094PC0012, \
    				 gts.1094PC0013 \
                 from variants" \
				 --header \
				 test.snpEff.vcf.db
    chrom	start	end	ref	alt	gene gts.1094PC0005	gts.1094PC0009	gts.1094PC0012	gts.1094PC0013
	chr1	30547	30548	T	G	FAM138A	./.	./.	./.	./.
    chr1	30859	30860	G	C	FAM138A	G/G	G/G	G/G	G/G
    chr1	30866	30869	CCT	C	FAM138A	CCT/CCT	CCT/CCT	CCT/C	CCT/CCT
    chr1	30894	30895	T	C	FAM138A	T/C	T/C	T/T	T/T
    chr1	30922	30923	G	T	FAM138A	./.	./.	./.	./.
    chr1	69269	69270	A	G	OR4F5	./.	./.	G/G	G/G
    chr1	69427	69428	T	G	OR4F5	T/T	T/T	T/T	T/T
    chr1	69510	69511	A	G	OR4F5	./.	./.	A/G	A/G
    chr1	69760	69761	A	T	OR4F5	A/A	A/T	A/A	A/A
    chr1	69870	69871	G	A	OR4F5	./.	G/G	G/G	G/G

Let's now get the genotype and the depth of aligned sequence observed for a
sample so that we can assess the confidence in the genotype::

	$ gemini query -q "select chrom, start, end, ref, alt, gene,
    				 gts.1094PC0005, \
					 gt_depths.1094PC0005, \
             	 from variants" test.snpEff.vcf.db
	chr1	30547	30548	T	G	FAM138A	./.	-1
	chr1	30859	30860	G	C	FAM138A	G/G	7
	chr1	30866	30869	CCT	C	FAM138A	CCT/CCT	8
	chr1	30894	30895	T	C	FAM138A	T/C	8
	chr1	30922	30923	G	T	FAM138A	./.	-1
	chr1	69269	69270	A	G	OR4F5	./.	-1
	chr1	69427	69428	T	G	OR4F5	T/T	2
	chr1	69510	69511	A	G	OR4F5	./.	-1
	chr1	69760	69761	A	T	OR4F5	A/A	1
	chr1	69870	69871	G	A	OR4F5	./.	-1

Now, we often want to focus only on variants where a given sample has a
specific genotype (e.g., looking for homozygous variants in family trios).
Unfortunately, we cannot directly do this in the SQL query, but the `gemini query`
tool has an option called `--gt-filter` that allows one to specify filters to
apply to the returned rows.  The rules followed in the `--gt-filter` option
follow Python syntax.  As an example, let's only return rows where sample
1094PC0012 is heterozygous.  In order to do this, we apply a filter to the
`gt_types` columns for this individual::

	$ gemini query -q "select chrom, start, end, ref, alt, gene,
    				 gts.1094PC0005, \
				 	 gts.1094PC0009, \
                 	 gts.1094PC0012, \
				 	 gts.1094PC0013 \
             	     from variants" \
			 	  --gt-filter "gt_types.1094PC0012 == HET" \
				  --header \
			 test.snpEff.vcf.db
	chrom	start	end	ref	alt	gene gts.1094PC0005	gts.1094PC0009	gts.1094PC0012	gts.1094PC0013
	chr1	30866	30869	CCT	C	FAM138A	CCT/CCT	CCT/CCT	CCT/C	CCT/CCT
	chr1	69510	69511	A	G	OR4F5	./.	./.	A/G	A/G

Now let's be a bit less restrictive and return variants where either sample
1094PC0012 is heterozygous or sample 1094PC0005 is homozygous for the reference
allele::

	$ gemini query -q "select chrom, start, end, ref, alt, gene,
					gts.1094PC0005, \
			 	 	gts.1094PC0009, \
             	 	gts.1094PC0012, \
			 	 	gts.1094PC0013 \
         	    	from variants" \
		 	  	 --gt-filter "gt_types.1094PC0012 == HET or \
				              gt_types.1094PC0005 == HOM_REF" \
			     --header \
		 test.snpEff.vcf.db
	chrom	start	end	ref	alt	gene gts.1094PC0005	gts.1094PC0009	gts.1094PC0012	gts.1094PC0013
	chr1	30859	30860	G	C	FAM138A	G/G	G/G	G/G	G/G
	chr1	30866	30869	CCT	C	FAM138A	CCT/CCT	CCT/CCT	CCT/C	CCT/CCT
	chr1	69427	69428	T	G	OR4F5	T/T	T/T	T/T	T/T
	chr1	69510	69511	A	G	OR4F5	./.	./.	A/G	A/G
	chr1	69760	69761	A	T	OR4F5	A/A	A/T	A/A	A/A

I changed my mind, let's restrict the above to those variants where sample
1094PC0012 must also be heterozygous::

	$ gemini query -q "select chrom, start, end, ref, alt, gene,
					gts.1094PC0005, \
		 	 		gts.1094PC0009, \
         	 		gts.1094PC0012, \
		 	 		gts.1094PC0013 \
     	    		from variants" \
	 	  	 	  --gt-filter "(gt_types.1094PC0012 == HET or \
			                   gt_types.1094PC0005 == HOM_REF) \
							   and \
							   (gt_types.1094PC0013 == HET)" \
		          --header \
	 test.snpEff.vcf.db
	 chrom	start	end	ref	alt	gene gts.1094PC0005	gts.1094PC0009	gts.1094PC0012	gts.1094PC0013
	 chr1	69510	69511	A	G	OR4F5	./.	./.	A/G	A/G


Adding your own, custom annotations to the ``GEMINI`` framework
===========================================================
It is inevitable that researchers will want to enhance the ``GEMINI``
framework with their own, custom annotations. ``GEMINI`` provides a
sub-command called ``annotate`` for exactly this purpose.  As long as
you provide a ``tabix``'ed annotation file in either BED or VCF format,
the ``annotate`` tool will, for each variant in the ``variants`` table,
screen for overlaps in your annotation file and update a new column in the
``variants`` table that you may specify on the command line.  This is best
illustrated by example.

Let's assume you have already created a ``GEMINI`` database of a VCF file
using the ``load`` module.

	$ gemini load -v my.vcf -t VEP my.db

Add a new column called "my_col" that tracks whether a given variant overlapped (1)
or did not overlap (0) intervals in your annotation file.

	# Add my_col to the database as a _boolean_
	$ gemini annotate -f my_annos.bed.gz -c my_col -t boolean my.db

	# Now query the results
	$ gemini query -q "select chrom, start, end, variant_id, my_col from variants" my.db | head -3
	chr22	16504479	16504480	1	1
	chr22	16504488	16504489	2	0
	chr22	16504490	16504491	3	1

Add a new column called "my_col" that counts the number of overlaps a given variant
has with intervals in your annotation file.

	# Add my_col to the database as a _count_
	$ gemini annotate -f my_annos.bed.gz -c my_col -t count my.db

	# Now query the results
	$ gemini query -q "select chrom, start, end, variant_id, my_col from variants" my.db | head -3
	chr22	16504479	16504480	1	2
	chr22	16504488	16504489	2	0
	chr22	16504490	16504491	3	1

Add a new column called "my_col" that creates a list of a specific column from the
annotation file for each given variant.

	# Add my_col to the database as a _list_, extracting the 4th
	# column from the annotation file to create the list
	$ gemini annotate -f my_annos.bed.gz -c my_col -t list -e 4 my.db

	# Now query the results
	$ gemini query -q "select chrom, start, end, variant_id, my_col from variants" my.db | head -3
	chr22	16504479	16504480	1	rs123,rs456
	chr22	16504488	16504489	2	None
	chr22	16504490	16504491	3	rs789


Extracting KEGG pathways for variants/genes/individuals
=======================================================

``GEMINI`` catalogs KEGG pathway information and using the
``pathways`` tool, one can extract pathway information for each
sample that has variants affecting a given gene.  The only requirement
is that we know what version of Ensembl genes were using by snpEff or VEP.
Currently, we expect version 66, 67, or 68.

    $ gemini pathways -v 66 chr22.low.exome.snpeff.100samples.vcf.db

We can also focus solely on loss-of-function mutations with the ``--lof``
argument:

    $ gemini pathways --lof -v 66 chr22.low.exome.snpeff.100samples.vcf.db | head
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


Filtering candidate LoF variants by transcript position or transcript type.
===========================================================================
Not all candidate LoF variants are created equal.  For example, a nonsense
(stop gain) variant impacting the first 5% of a polypeptide is far more
likely to be deleterious than one affecting the last 5%. (For an empirical
analysis of this in the human genome, see Fig 1C in [MacArthur et al, 2012](http://www.sciencemag.org/content/335/6070/823)). Assuming you've
annotated your VCF with snpEff v3.0+, the ``lof_sieve`` tool
reports the fractional position (e.g. 0.05 for the first 5%) of the mutation
in the amino acid sequence. In addition, it also reports the predicted
function of the transcript so that one can segregate candidate LoF
variants that affect protein_coding transcripts from processed RNA, etc.

    $ gemini lof_sieve chr22.low.exome.snpeff.100samples.vcf.db
    chrom	start	end	ref	alt	highest_impact	aa_change	var_trans_pos	trans_aa_length	var_trans_pct	sample	genotype	gene	transcript	trans_type
    chr22	17072346	17072347	C	T	stop_gain	W365*	365	557	0.655296229803	NA19327	C|T	CCT8L2	ENST00000359963	protein_coding
    chr22	17072346	17072347	C	T	stop_gain	W365*	365	557	0.655296229803	NA19375	T|C	CCT8L2	ENST00000359963	protein_coding
    chr22	17072346	17072347	C	T	stop_gain	W365*	365	557	0.655296229803	NA19431	T|C	CCT8L2	ENST00000359963	protein_coding
    chr22	17129539	17129540	C	T	splice_donor	None	None	None	None	NA18964	T|C	TPTEP1	ENST00000383140	lincRNA
    chr22	17129539	17129540	C	T	splice_donor	None	None	None	None	NA19675	T|C	TPTEP1	ENST00000383140	lincRNA
    chr22	17140745	17140746	A	G	splice_donor	None	None	None	None	NA19223	A|G	ANKRD62P1	ENST00000456726	lincRNA


Identifying potential compound heterozygotes
===========================================================================
Many recessive disorders are caused by compound heterozygotes. Unlike canonical
recessive sites where the same recessive allele is inherited from both parents
at the _same_ site in the gene, compund heterozygotes occur when
the individual's phenotype is caused by two heterogeneous recessive alleles at
_different_ sites in a particular gene.

So basically, we are looking for two (typically loss-of-function (LoF))
heterozygous variants impacting the same gene at different loci.  The
complicating factor is that this is _recessive_ and as such, we must also
require that the consequential alleles at each heterozygous site were
inherited on different chromosomes (one from each parent).  As such, in order
to use this tool, we require that all variants are phased.  Once this has been
done, the `comp_hets` tool will provide a report of candidate compund
heterozygotes for each sample/gene.

For example:

	$ gemini comp_hets chr22.low.exome.snpeff.100samples.vcf.db
	sample	gene	het1	het2
	NA19002	GTSE1	chr22,46722400,46722401,G,A,G|A,stop_gain,exon_22,0.005,1	chr22,46704499,46704500,C,A,A|C,stop_gain,exon_22,0.005,0

This indicates that sample NA19002 has a candidate compund heterozygote in
GTSE1.  The two hets are reported using the following structure:
(chrom,start,end,ref,alt,genotype,impact,exon,AAF,in_dbsnp).

By default, all coding variants are explored.  However, one may want to
restrict the analysis to LoF variants.

	$ gemini comp_hets --only_lof chr22.low.exome.snpeff.100samples.vcf.db

Extracting variants from specific regions or genes
==================================================

``GEMINI`` allows one to extract variants that fall within specific genomic coordinates as follows

	$ gemini region --reg chr1:100-200 my.db

Or, one can extract variants based on a specific gene name.

	$ gemini region --gene PTPN22 my.db


Conducting analyses on genome "windows".
========================================

``GEMINI`` includes a convenient tool for computing variation metrics across genomic windows (fixed and sliding).
Here are a few examples to whet your appetite.  If you're still hungry, email us.

Compute the average nucleotide diversity for all variants found in non-overlapping, 50Kb windows.

	$ gemini windower -w 50000 -s 0 -t nucl_div -o mean my.db

Compute the average nucleotide diversity for all variants found in 50Kb windows that overlap by 10kb.

	$ gemini windower -w 50000 -s 10000 -t nucl_div -o mean my.db


I forget the database table and column names.  Help!
====================================================

Just use the `db_info` tool.  It will report all of the tables and all of the
columns / types in each table:

	$ gemini db_info test.db
	table_name          column_name                   type
	variants            chrom                         text
	variants            start                         integer
	variants            end                           integer
	variants            variant_id                    integer
	variants            anno_id                       integer
	variants            ref                           text
	variants            alt                           text
	variants            qual                          float
	variants            filter                        text
	variants            type                          text
	variants            sub_type                      text
	variants            gts                           blob
	variants            gt_types                      blob
	variants            gt_phases                     blob
	variants            gt_depths                     blob
	variants            call_rate                     float
	variants            in_dbsnp                      bool
	variants            rs_ids                        text
	variants            in_omim                       bool
	variants            clin_sigs                     text
	variants            cyto_band                     text
	variants            rmsk                          text
	variants            in_cpg_island                 bool
	variants            in_segdup                     bool
	variants            is_conserved                  bool
	variants            num_hom_ref                   integer
	variants            num_het                       integer
	variants            num_hom_alt                   integer
	variants            num_unknown                   integer
	variants            aaf                           float
	variants            hwe                           float
	variants            inbreeding_coeff              float
	variants            pi                            float
	variants            recomb_rate                   float
	variants            gene                          text
	variants            transcript                    text
	variants            is_exonic                     bool
	variants            is_coding                     bool
	variants            is_lof                        bool
	variants            exon                          text
	variants            codon_change                  text
	variants            aa_change                     text
	variants            aa_length                     text
	variants            biotype                       text
	variants            impact                        text
	variants            impact_severity               text
	variants            polyphen_pred                 text
	variants            polyphen_score                float
	variants            sift_pred                     text
	variants            sift_score                    float
	variants            anc_allele                    text
	variants            rms_bq                        float
	variants            cigar                         text
	variants            depth                         integer
	variants            strand_bias                   float
	variants            rms_map_qual                  float
	variants            in_hom_run                    integer
	variants            num_mapq_zero                 integer
	variants            num_alleles                   integer
	variants            num_reads_w_dels              float
	variants            haplotype_score               float
	variants            qual_depth                    float
	variants            allele_count                  integer
	variants            allele_bal                    float
	variants            in_hm2                        bool
	variants            in_hm3                        bool
	variants            is_somatic
	variants            in_esp                        bool
	variants            aaf_esp_ea                    float
	variants            aaf_esp_aa                    float
	variants            aaf_esp_all                   float
	variants            exome_chip                    bool
	variants            in_1kg                        bool
	variants            aaf_1kg_amr                   float
	variants            aaf_1kg_asn                   float
	variants            aaf_1kg_afr                   float
	variants            aaf_1kg_eur                   float
	variants            aaf_1kg_all                   float
	variants            grc                           text
	variants            gms_illumina                  float
	variants            gms_solid                     float
	variants            gms_iontorrent                float
	variants            encode_tfbs
	variants            encode_consensus_gm12878      text
	variants            encode_consensus_h1hesc       text
	variants            encode_consensus_helas3       text
	variants            encode_consensus_hepg2        text
	variants            encode_consensus_huvec        text
	variants            encode_consensus_k562         text
	variants            encode_segway_gm12878         text
	variants            encode_segway_h1hesc          text
	variants            encode_segway_helas3          text
	variants            encode_segway_hepg2           text
	variants            encode_segway_huvec           text
	variants            encode_segway_k562            text
	variants            encode_chromhmm_gm12878       text
	variants            encode_chromhmm_h1hesc        text
	variants            encode_chromhmm_helas3        text
	variants            encode_chromhmm_hepg2         text
	variants            encode_chromhmm_huvec         text
	variants            encode_chromhmm_k562          text
	variant_impacts     variant_id                    integer
	variant_impacts     anno_id                       integer
	variant_impacts     gene                          text
	variant_impacts     transcript                    text
	variant_impacts     is_exonic                     bool
	variant_impacts     is_coding                     bool
	variant_impacts     is_lof                        bool
	variant_impacts     exon                          text
	variant_impacts     codon_change                  text
	variant_impacts     aa_change                     text
	variant_impacts     aa_length                     text
	variant_impacts     biotype                       text
	variant_impacts     impact                        text
	variant_impacts     impact_severity               text
	variant_impacts     polyphen_pred                 text
	variant_impacts     polyphen_score                float
	variant_impacts     sift_pred                     text
	variant_impacts     sift_score                    float
	samples             sample_id                     integer
	samples             name                          text
	samples             family_id                     integer
	samples             paternal_id                   integer
	samples             maternal_id                   integer
	samples             sex                           text
	samples             phenotype                     text
	samples             ethnicity                     text
