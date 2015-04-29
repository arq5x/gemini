###############################
Querying the GEMINI database
###############################


The real power in the ``GEMINI`` framework lies in the fact that all of your
genetic variants have been stored in a convenient database in the context of a
wealth of genome annotations that facilitate variant interpretation.  The
expressive power of SQL allows one to pose intricate questions of one's variation
data.

.. note::

    If you are unfamiliar with SQL, `sqlzoo <http://sqlzoo.net/>`_ has a decent
    online tutorial describing the basics.  Really all you need to learn is the
    SELECT statement, and the examples below will give you a flavor of how to
    compose base SQL queries against the GEMINI framework.


==============================
Basic queries
==============================

GEMINI has a specific tool for querying a gemini database that has been loaded
using the ``gemini load`` command.  That's right, the tool is called
``gemini query``. Below are a few basic queries that give you a sense of how to
interact with the gemini database using the ``query`` tool.

1. Extract all transitions with a call rate > 95%

.. code-block:: bash

    $ gemini query -q "select * from variants \
                          where sub_type = 'ts' \
                          and call_rate >= 0.95" my.db

2. Extract all loss-of-function variants with an alternate allele frequency < 1%:

.. code-block:: bash

    $ gemini query -q "select * from variants \
                          where is_lof = 1 \
                          and aaf >= 0.01" my.db

3. Extract the nucleotide diversity for each variant:

.. code-block:: bash

    $ gemini query -q "select chrom, start, end, pi from variants" my.db


4. Combine ``GEMINI`` with ``bedtools`` to compute nucleotide diversity estimates across 100kb windows:

.. code-block:: bash

    $ gemini query -q "select chrom, start, end, pi from variants \
                          order by chrom, start, end" my.db | \
      bedtools map -a hg19.windows.bed -b - -c 4 -o mean


===========================================================
Selecting sample genotypes
===========================================================

The above examples illustrate *ad hoc* queries that do not request or filter
upon the genotypes of individual samples.  Since ``GEMINI`` stores the genotype
information for each variant in compressed arrays that are stored as BLOBs
in the database, standard SQL queries cannot directly access individual
genotypes. However, we have enhanced the SQL syntax to support such queries
with C "struct-like" access.  For example, to retrieve the alleles for a given
sample's (in this case, sample 1094PC0009), one would add ``gts.1094PC0009``
to the select statement.

Here is an example of selecting the genotype alleles for four
different samples (note the examples below use the test.snpEff.vcf.db
file that is created in the ./test directory when you run the
`bash master-test.sh` command as described above):

.. code-block:: bash

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


You can also add a header so that you can keep track of who's who:

.. code-block:: bash

    $ gemini query -q "select chrom, start, end, ref, alt, gene, \
                              gts.1094PC0005, \
                              gts.1094PC0009, \
                              gts.1094PC0012, \
                              gts.1094PC0013 \
                       from variants" \
                       --header test.snpEff.vcf.db

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
sample so that we can assess the confidence in the genotype:

.. code-block:: bash

	$ gemini query -q "select chrom, start, end, ref, alt, gene,
                              gts.1094PC0005, \
                              gt_depths.1094PC0005 \
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

===========================================================
Selecting sample genotypes based on "wildcards".
===========================================================
The above examples demonstrate how one can select individual sample genotype
information by explicitly listing each column and sample that one wishes to see.
Obviously, this can become tedious when a project involves hundreds or thousands of samples 
--- if you wanted to see genotype information for the 345 of 1145 affected samples in your study,
you would have to type each and every column.sample name out. Brutal.

-----------------------------------------------------------
The "*" wildcard
-----------------------------------------------------------

To get around this, one can bulk-select sample genotype information using "wildcards". The column and the wildcard must each be surrounded with parentheses and separated by a period. The "*" is a shortcut (wildcard) meaning "all samples".

.. note::

  The syntax for SELECTing genotype columns using a wildcard is ``(COLUMN).(WILDCARD)``.

For example, a shortcut to reporting the genotype for *all* samples (in this case 4) in the study, one could do the following:

.. code-block:: bash

  $ gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(*) \
                              from variants" extended_ped.db
  chrom start end ref alt gene  gts.M10475  gts.M10478  gts.M10500  gts.M128215
  chr10 1142207 1142208 T C WDR37 C/C C/C C/C C/C
  chr10 48003991  48003992  C T ASAH2C  T/T C/T C/T C/C
  chr10 52004314  52004315  T C ASAH2 ./. ./. C/C C/C
  chr10 52497528  52497529  G C ASAH2B  ./. C/C C/C ./.
  chr10 126678091 126678092 G A CTBP2 G/G G/G G/G G/A
  chr10 135210790 135210791 T C MTG1.1  T/T C/C C/C T/T
  chr10 135336655 135336656 G A SPRN  ./. A/A ./. A/A
  chr10 135369531 135369532 T C SYCE1 T/T T/C T/C T/T
  chr16 72057434  72057435  C T DHODH C/T C/C C/C C/C

-----------------------------------------------------------
Wildcards based on sample attributes
-----------------------------------------------------------

To report the genotypes for solely those samples that are affected (phenotype == 2) with the phenotype in question, one could do the following:

.. code-block:: bash

  $ gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(phenotype==2) \
                              from variants" extended_ped.db
  chrom start end ref alt gene  gts.M10478  gts.M10500
  chr10 1142207 1142208 T C WDR37 C/C C/C
  chr10 48003991  48003992  C T ASAH2C  C/T C/T
  chr10 52004314  52004315  T C ASAH2 ./. C/C
  chr10 52497528  52497529  G C ASAH2B  C/C C/C
  chr10 126678091 126678092 G A CTBP2 G/G G/G
  chr10 135210790 135210791 T C MTG1.1  C/C C/C
  chr10 135336655 135336656 G A SPRN  A/A ./.
  chr10 135369531 135369532 T C SYCE1 T/C T/C
  chr16 72057434  72057435  C T DHODH C/C C/C

One can add multiple wildcard criteria as well:

.. code-block:: bash

  $ gemini query --header -q "select chrom, start, end, ref, alt, gene, (gts).(phenotype==1 and hair_color=='blue') \
                              from variants" extended_ped.db
  chrom start end ref alt gene  gts.M128215
  chr10 1142207 1142208 T C WDR37 C/C
  chr10 48003991  48003992  C T ASAH2C  C/C
  chr10 52004314  52004315  T C ASAH2 C/C
  chr10 52497528  52497529  G C ASAH2B  ./.
  chr10 126678091 126678092 G A CTBP2 G/A
  chr10 135210790 135210791 T C MTG1.1  T/T
  chr10 135336655 135336656 G A SPRN  A/A
  chr10 135369531 135369532 T C SYCE1 T/T
  chr16 72057434  72057435  C T DHODH C/C

===========================================================
``--gt-filter`` Filtering on genotypes
===========================================================

Now, we often want to focus only on variants where a given sample has a
specific genotype (e.g., looking for homozygous variants in family trios).
Unfortunately, we cannot directly do this in the SQL query, but the `gemini query`
tool has an option called `--gt-filter` that allows one to specify filters to
apply to the returned rows.  The rules followed in the `--gt-filter` option
follow Python syntax.


.. tip::

    As you will see from the examples below, appropriate use of the --gt-filter
    option will allow you to compose queries that return variants meeting
    inheritance patterns that are relevant to the disease model of interest
    in your study.

As an example, let's only return rows where sample
1094PC0012 is heterozygous.  In order to do this, we apply a filter to the
`gt_types` columns for this individual:

.. code-block:: bash

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
allele:

.. code-block:: bash

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


Eh, I changed my mind, let's restrict the above to those variants where sample
1094PC0012 must also be heterozygous:

.. code-block:: bash

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

===========================================================
``--gt-filter`` Wildcard filtering on genotype columns.
===========================================================

Many times, we want to be able to apply the same rule to multiple samples
without having to enter the rule over and over again for each sample. For example,
let's imaging there are 100 samples in your study and you only want to report variants
where every sample has an observed alignment depth of at least 20 reads. Traditionally,
one would have enter each of the 100 samples from the command line as follows:

.. code-block:: bash

  $ gemini query -q "select chrom, start, end, ref, alt, gene from variants" \
                       --gt-filter "gt_depths.sample1 >= 20 and \
                                    gt_depths.sample2 >= 20 and \
                                    gt_depths.sample3 >= 20 and \
                                    ...
                                    gt_depths.sample100 >= 20" \
                       extended_ped.db


-----------------------------------------------------------
The basics.
-----------------------------------------------------------

Obviously, this is deeply painful. Thankfully, there an option allowing the use of "wildcards" to prevent this. 

.. note::

  The syntax of the wildcard ``--gt-filters`` is ``(COLUMN).(SAMPLE_WILDCARD).(SAMPLE_WILDCARD_RULE).(RULE_ENFORCEMENT)``.

For example, using wildcards, the above could be converted to:

.. code-block:: bash

  $ gemini query -q "select chrom, start, end, ref, alt, gene from variants" \
                       --gt-filter "(gt_depths).(*).(>=20).(all)" \
                       extended_ped.db


Obviously, this makes things much simpler.


-----------------------------------------------------------
The ``all`` operator
-----------------------------------------------------------

One can also apply wildcards that select samples based on the values in specific 
columns in the ``samples`` table. For example, let's imagine we wanted to require that variants
are returned only in cases where *ALL* the affected individuals in the study (i.e., the ``phenotype`` column in the ``samples`` table is equal to ``2``) have non-reference genotypes. We could do the following:

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, (gts).(phenotype==2) \
                     FROM variants" \
                 --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all)" \
                 extended_ped.db
  chrom start end ref alt gene  gts.M10478  gts.M10500
  chr10 1142207 1142208 T C WDR37 C/C C/C
  chr10 48003991  48003992  C T ASAH2C  C/T C/T
  chr10 52004314  52004315  T C ASAH2 ./. C/C
  chr10 52497528  52497529  G C ASAH2B  C/C C/C
  chr10 135210790 135210791 T C MTG1.1  C/C C/C
  chr10 135336655 135336656 G A SPRN  A/A ./.
  chr10 135369531 135369532 T C SYCE1 T/C T/C

Or perhaps we wanted to be more restrictive. We could also enforce that the affected individuals also had at least 20 aligned reads at such variant sites:

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, \
                            (gts).(phenotype==2), (gt_depths).(phenotype==2) \
                     FROM variants" \
                 --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) \
                              and \
                              (gt_depths).(phenotype==2).(>=20).(all)" \
                 extended_ped.db
  chrom start end ref alt gene  gts.M10478  gts.M10500  gt_depths.M10478  gt_depths.M10500
  chr10 1142207 1142208 T C WDR37 C/C C/C 29  24
  chr10 48003991  48003992  C T ASAH2C  C/T C/T 38  44
  chr10 135369531 135369532 T C SYCE1 T/C T/C 32  21

-----------------------------------------------------------
The ``any`` operator
-----------------------------------------------------------

The examples provided thus far have deomnstrated how to enforce that *ALL* of the samples meeting specific criteria meet the same genotype column restrictions. However, clearly there are cases where one would want to be less restrictive. In order to accomondate such queries, there are three other enforcement operators allowed: ``any``, ``none``, and ``count``.

For example, perhaps we want to relax the above query a bit and only require that at least one (i.e., ``any``) of the affected samples has depth > 20:

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, \
                            (gts).(phenotype==2), (gt_depths).(phenotype==2) \
                     FROM variants" \
                 --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) \
                              and \
                              (gt_depths).(phenotype==2).(>=20).(any)" \
                 extended_ped.db

-----------------------------------------------------------
The ``none`` operator
-----------------------------------------------------------

Or we enforce that ``none`` of the affected samples have depth less than 10:

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, \
                            (gts).(phenotype==2), (gt_depths).(phenotype==2) \
                     FROM variants" \
                 --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) \
                              and \
                              (gt_depths).(phenotype==2).(<10).(none)" \
                 extended_ped.db

-----------------------------------------------------------
The ``count`` operator
-----------------------------------------------------------

Finally, we could enforce that at most 2 of all the samples in the study have depths < 10:

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, \
                            (gts).(phenotype==2), (gt_depths).(phenotype==2) \
                     FROM variants" \
                 --gt-filter "(gt_types).(phenotype==2).(!=HOM_REF).(all) \
                              and \
                              (gt_depths).(phenotype==2).(<10).(count <= 2)" \
                 extended_ped.db

.. note::

  The ``count`` operator allows the following comparisons: ``>``, ``>=``, ``<``, ``<=``, ``==``, ``<>``, and ``!=``.

-----------------------------------------------------------
Example scenario.
-----------------------------------------------------------
One usage of the wildcard functionality is screening for variants
that are present in affected individuals yet absent from unaffected individuals (this is obviously an idealized scenario).

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, (gts).(*) \
                     FROM variants" \
                 --gt-filter "(gt_types).(phenotype==1).(==HOM_REF).(all) \
                              and \
                              (gt_depths).(phenotype==2).(==HOM_REF).(none)" \
                 extended_ped.db


-----------------------------------------------------------
More complexity.
-----------------------------------------------------------
One can also combine wildcard filters with "basic" genotype filters
to create more complex logic.

.. code-block:: bash

  $ gemini query --header \
                 -q "SELECT chrom, start, end, ref, alt, gene, (gts).(*) \
                     FROM variants" \
                 --gt-filter "((gt_types).(phenotype==1).(==HOM_REF).(all) \
                                and \
                               (gt_depths).(phenotype==2).(==HOM_REF).(none)) \
                              or gt.NA12878 == HET" \
                 extended_ped.db


-----------------------------------------------------------
Other details.
-----------------------------------------------------------

The system is fairly flexible in that it allows one to wildcard-select samples based on custom columns 
that have been added to the ``samples`` table based upon a custom PED file.  For example, let's imaging our
custom PED file had an extra column defining the hair color of each sample. We coukd use that to restrict interesting variants to those where samples with blue hair were heterozygous:

.. code-block:: bash

  $ gemini query -q "SELECT chrom, start, end, ref, alt, gene FROM variants" \
                 --gt-filter "(gt_types).(hair_color='blue').(==HET).(all)" \
                 extended_ped.db

Or possibly, you want to stratify based on sub-population:

.. code-block:: bash

  $ gemini query -q "SELECT chrom, start, end, ref, alt, gene FROM variants" \
                 --gt-filter "(gt_types).(population='CEU').(==HET).(all) \
                               and \
                              (gt_types).(population='YRI').(==HOM_ALT).(any)" \
                 extended_ped.db


One can also base the wildcard on multiple criteria (in this case, brown hair and affected):

.. code-block:: bash

  $ gemini query -q "SELECT chrom, start, end, ref, alt, gene FROM variants" \
                 --gt-filter "(gt_types).(hair_color=='brown' and phenotype==2).(!= HET).(all)" \
                 extended_ped.db

Lastly, wildcards can, of course, be combined with non-wildcard criteria:

.. code-block:: bash

  $ gemini query -q "SELECT chrom, start, end, gene FROM variants" \
                 --gt-filter "(gt_types).(hair_color=='brown' and phenotype==2).(!= HET).(all) \
                               and \
                              gt_types.M128215 == HOM_REF" \
                 extended_ped.db

Hopefully this gives you a sense of what you can do with the "wildcard" genotype filter functionality.


=============================================================
``--show-samples`` Finding out which samples have a variant
=============================================================
While exploring your data you might hit on a set of interesting variants and want to know
which of your samples have that variant in them. You can display the samples containing
a variant with the --show-sample-variants flag:

.. code-block:: bash

   	$ gemini query --header --show-samples -q "select chrom, start, end, ref, alt \
                                        from variants where is_lof=1 limit 5" test.query.db

	chrom	start	end	ref	alt	variant_samples	het_samples	hom_alt_samples
	chr1	874815	874816	C	CT	1478PC0006B,1478PC0007B,1478PC0010,1478PC0013B,1478PC0022B,1478PC0023B,1478PC0025,1719PC0007,1719PC0009,1719PC0010,1719PC0022	1478PC0006B,1478PC0007B,1478PC0010,1478PC0013B,1478PC0022B,1478PC0023B,1719PC0007,1719PC0009,1719PC0010	1478PC0025,1719PC0022
	chr1	1140811	1140813	TC	T	1478PC0011	1478PC0011
	chr1	1219381	1219382	C	G	1719PC0012	1719PC0012
	chr1	1221487	1221490	CAA	C	1478PC0004	1478PC0004

variant_samples is a list of all of the samples with a variant, het_samples is the subset
of those heterozygous for the variant and hom_alt_samples is the subset homozygous for
the variant.

============================================================================
``--show-samples --format sampledetail`` Provide a flattened view of samples
============================================================================
If you'd like to be able to export variants plus associated sample metadata into
a downstream tool like R or pandas for additional exploration, adding the
``--format sampledetail`` command flattens all found samples and attached
metadata information:

.. code-block:: bash

    $ gemini query --header --format sampledetail --show-samples \
       -q "select chrom, start, ref from variants where is_lof=1 limit 1" test.query.db

    chrom	start	ref	family_id	name	paternal_id	maternal_id	sex	phenotype
    chr1	30547	T	0	1478PC0016	0	0	-9	-9
    chr1	30547	T	0	1719PC0007	0	0	-9	-9
    chr1	30547	T	0	1719PC0009	0	0	-9	-9" > exp

The denormalized results contain a row for each sample associated with a
variant, along with information from the sample table. This provides a way to
join sample information with a variant query.

=============================================================
``--show-families`` Finding out which families have a variant
=============================================================
This works exactly like ``--show-samples`` except lists all of the families with a
variant instead of the individual samples.

===================================================
``--region`` Restrict a query to a specified region
===================================================
If you are only interested in a specific region, you can restrict queries to
that region using the ``--region`` tool.

.. code-block:: bash

   $ gemini query --region chr1:30859-30900 -q "select chrom, start, end, ref, alt \
                from variants"  test1.snpeff.db
   chr1	30859	30860	G	C

=========================================================
``--sample-filter`` Restrict a query to specified samples
=========================================================
The ``--sample-filter`` option allows you to select samples that a variant
must be in by doing a SQL query on the samples table. For example if you
wanted to show the set of variants that appear in all samples with
a phenotype status of 2, you could do that query with:

.. code-block:: bash

   $ gemini query --sample-filter "phenotype=2" -q "select gts, gt_types from variants" test.family.db
   T/T,T/T,T/C,T/T,T/T,T/T,T/T,T/T,C/C	0,0,1,0,0,0,0,0,3	1_kid,3_kid	1_kid	3_kid
   T/T,T/T,T/C,T/T,T/T,T/C,T/T,T/T,T/C	0,0,1,0,0,1,0,0,1	1_kid,2_kid,3_kid	1_kid,2_kid,3_kid
   T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/C	0,0,0,0,0,0,0,0,1	3_kid	3_kid

By default --sample-filter will show the variant if at least one sample contains the
variant. You can change this behavior by using the ``--in`` option along with
``--sample-filter``. ``--in all`` will return a variant if all samples matching
the query have the variant. ``in none`` will return a variant if the variant
does not appear in any of the matching samples. ``--in only`` will return a variant
if the variant is only in the matching samples and not in any of the non-matching
samples. ``--in only all`` will show all of the variant which are in all of the
matching samples and not in any of the non-matching samples.

The ``--family-wise`` flag applies the ``--sample-filter`` and ``--in`` behavior
on a family-wise basis. For example to show all variants that are only in samples
with a phenotype status of 2 in at least one family:

.. code-block:: bash

   $ gemini query --family-wise --in only all --sample-filter "phenotype=2" -q "select gts, gt_types from variants" test.family.db
   T/T,T/T,T/C,T/T,T/T,T/T,T/T,T/T,C/C	0,0,1,0,0,0,0,0,3	1_kid,3_kid	1_kid	3_kid
   T/T,T/T,T/C,T/T,T/T,T/C,T/T,T/T,T/C	0,0,1,0,0,1,0,0,1	1_kid,2_kid,3_kid	1_kid,2_kid,3_kid
   T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/C	0,0,0,0,0,0,0,0,1	3_kid	3_kid


You can also specify that a variant passes this filter in multiple families with
the ``--min-kindreds`` option. So if you want to do the same query above, but restrict it
such that at least three families have to pass the filter:

.. code-block:: bash

   $ gemini query --min-kindreds 3 --family-wise --in only all --sample-filter "phenotype=2" -q "select gts, gt_types from variants" test.family.db
   T/T,T/T,T/C,T/T,T/T,T/C,T/T,T/T,T/C	0,0,1,0,0,1,0,0,1	1_kid,2_kid,3_kid	1_kid,2_kid,3_kid


If the PED file you loaded has extra fields in it, those will also work with the
``--sample-filter`` option. For example if you had a ``hair_color`` extended field,
you could query on that as well as phenotype:

.. code-block:: bash

   $ gemini query  --in only all --sample-filter "phenotype=1 and hair_color='blue'" -q "select gts, gt_types from variants" extended_ped.db
   G/G,G/G,G/G,G/A	0,0,0,1	M128215	M128215


=============================================================
``--sample-delim`` Changing the sample list delimiter
=============================================================
One can modify the default comma delimiter used by the ``--show-samples``
option through the use of the ``--sample-delim`` option.  For example, to use
a semi-colon instead of a comma, one would do the following:

.. code-block:: bash

    $ gemini query --header --show-samples --sample-delim ";" \
                   -q "select chrom, start, end, ref, alt \
                       from variants where is_lof=1 limit 5" test.query.db

  chrom start end ref alt variant_samples het_samples hom_alt_samples
  chr1  874815  874816  C CT  1478PC0006B;1478PC0007B;1478PC0010,1478PC0013B;1478PC0022B;1478PC0023B;1478PC0025;1719PC0007;1719PC0009;1719PC0010;1719PC0022 1478PC0006B;1478PC0007B;1478PC0010;1478PC0013B;1478PC0022B;1478PC0023B;1719PC0007;1719PC0009;1719PC0010 1478PC0025;1719PC0022
  chr1  1140811 1140813 TC  T 1478PC0011  1478PC0011
  chr1  1219381 1219382 C G 1719PC0012  1719PC0012
  chr1  1221487 1221490 CAA C 1478PC0004  1478PC0004


===========================================================
``--format`` Reporting query output in an alternate format.
===========================================================
The results of GEMINI queries can automatically be formatted for use with
other programs using the --format command. Supported alternative
formats are JSON and TPED (Transposed PED) format.

Reporting query output in JSON format may enable
HTML/Javascript apps to query GEMINI and retrieve
the output in a format that is amenable to web development protocols.

Here is a basic query:

.. code-block:: bash

  $ gemini query -q "select chrom, start, end from variants" my.db | head
  chr1  10067 10069
  chr1  10230 10231
  chr1  12782 12783
  chr1  13109 13110
  chr1  13115 13116
  chr1  13117 13118
  chr1  13272 13273
  chr1  13301 13302
  chr1  13416 13417
  chr1  13417 13418

To report in JSON format, use the ``--format json`` option. For example:


.. code-block:: bash

  $ gemini query --format json -q "select chrom, start, end from variants" my.db | head
  {"chrom": "chr1", "start": 10067, "end": 10069}
  {"chrom": "chr1", "start": 10230, "end": 10231}
  {"chrom": "chr1", "start": 12782, "end": 12783}
  {"chrom": "chr1", "start": 13109, "end": 13110}
  {"chrom": "chr1", "start": 13115, "end": 13116}
  {"chrom": "chr1", "start": 13117, "end": 13118}
  {"chrom": "chr1", "start": 13272, "end": 13273}
  {"chrom": "chr1", "start": 13301, "end": 13302}
  {"chrom": "chr1", "start": 13416, "end": 13417}
  {"chrom": "chr1", "start": 13417, "end": 13418}


If you would to use tools such as PLINK that use
the PED format, you can dump out a set of variants
matching any query in TPED (Transposed PED) format
by adding the ``--tped`` flag to your query:

.. code-block:: bash

    $ gemini query --format tped -q "select * from variants where chrom=10" test4.snpeff.db
    10 rs10794716 0 1142207 C/C C/C C/C C/C
    10 rs142685947 0 48003991 T/T C/T C/T C/C
    10 rs2842123 0 52004314 ./. ./. C/C C/C
    10 rs4935178 0 52497528 ./. C/C C/C ./.
    16 rs201947120 0 72057434 C/T C/C C/C C/C
    10 rs73373169 0 126678091 G/G G/G G/G G/A
    10 rs2265637 0 135210790 T/T C/C C/C T/T
    10 rs6537611 0 135336655 ./. A/A ./. A/A
    10 rs3747881 0 135369531 T/T T/C T/C T/T

You can pass --header to get a header to see which samples have which
variant. To use the TPED format you also need to generate a corresponing TFAM
file from your data as well, which you can get from the GEMINI dump tool:

.. code-block:: bash

    $ gemini dump  --tfam test4.snpeff.db > obs
    None    M10475  None    None    None    None
    None    M10478  None    None    None    None
    None    M10500  None    None    None    None
    None    M128215 None    None    None    None

===========================================================
``--carrier-summary-by-phenotype`` Summarize carrier status
===========================================================
For prioritizing variants sometimes it is useful to have summary counts of
the carrier status for all samples with a variant stratified across a phenotype.
``--carrier-summary-by-phenotype`` takes a column in the samples table that you
want to summarize the carrier status of and adds a set of counts of
carrier/non-carrier status for each phenotype in the given column. For example,
to get a summary of how a set of variants segregate with affected status:

.. code-block:: bash

	$ gemini query --show-samples --carrier-summary-by-phenotype affected --header -q "select chrom, start, ref, alt, gt_types from variants" extended_ped_test.db
	chrom	start	ref	alt	gt_types	variant_samples	het_samples	hom_alt_samples	unaffected_carrier	affected_carrier	unaffected_noncarrier	affected_noncarrier	unknown
	chr10	1142207	T	C	3,3,3,3	M10475,M10478,M10500,M128215		M10475,M10478,M10500,M128215	2	2	0	0	0
	chr10	48003991	C	T	3,1,1,0	M10475,M10478,M10500	M10478,M10500	M10475	1	2	1	0	0
	chr10	52004314	T	C	2,2,3,3	M10500,M128215		M10500,M128215	1	1	0	0	2
	chr10	52497528	G	C	2,3,3,2	M10478,M10500		M10478,M10500	0	2	0	0	2
	chr16	72057434	C	T	1,0,0,0	M10475	M10475		1	0	1	2	0
	chr10	126678091	G	A	0,0,0,1	M128215	M128215		1	0	1	2	0
	chr10	135210790	T	C	0,3,3,0	M10478,M10500		M10478,M10500	0	2	2	0	0
	chr10	135336655	G	A	2,3,2,3	M10478,M128215		M10478,M128215	1	1	0	0	2
	chr10	135369531	T	C	0,1,1,0	M10478,M10500	M10478,M10500		0	2	2	0	0

Or if you have another phenotypic feature you are interested in summarizing,
like hair color:

.. code-block:: bash

	$ gemini query --show-samples --carrier-summary-by-phenotype hair_color --header -q "select chrom, start, ref, alt, gt_types from variants" extended_ped.db
	chrom	start	ref	alt	gt_types	variant_samples	het_samples	hom_alt_samples	blue_carrier	brown_carrier	purple_carrier	blue_noncarrier	brown_noncarrier	purple_noncarrier	unknown
	chr10	1142207	T	C	3,3,3,3	M10475,M10478,M10500,M128215		M10475,M10478,M10500,M128215	1	2	1	0	0	0	0
	chr10	48003991	C	T	3,1,1,0	M10475,M10478,M10500	M10478,M10500	M10475	0	2	1	1	0	0	0
	chr10	52004314	T	C	2,2,3,3	M10500,M128215		M10500,M128215	1	0	1	0	0	0	2
	chr10	52497528	G	C	2,3,3,2	M10478,M10500		M10478,M10500	0	1	1	0	0	0	2
	chr16	72057434	C	T	1,0,0,0	M10475	M10475		0	1	0	1	1	1	0
	chr10	126678091	G	A	0,0,0,1	M128215	M128215		1	0	0	0	2	1	0
	chr10	135210790	T	C	0,3,3,0	M10478,M10500		M10478,M10500	0	1	1	1	1	0	0
	chr10	135336655	G	A	2,3,2,3	M10478,M128215		M10478,M128215	1	1	0	0	0	0	2
	chr10	135369531	T	C	0,1,1,0	M10478,M10500	M10478,M10500		0	1	1	1	1	0	0

#################################
Querying the gene tables
#################################
The gene tables viz. ``gene_detailed table`` and the ``gene_summary table`` have been built on version 73 of the ensembl genes. The column specifications are
available at :doc:`database_schema`. These tables contain gene specific information e.g. gene synonyms, RVIS percentile scores(Petrovski et.al 2013), strand specifications, cancer gene census etc. While the former is more detailed, the later lacks transcript wise information and summarizes some aspects of the former. For e.g. while the gene_detailed table lists all transcripts of a gene with their start and end co-ordinates, the gene_summary table reports only the minimum start and maximum end co-ordinates of the gene transcripts. The ``chrom``, ``gene`` and the ``transcript`` columns of the gene tables may be used to join on the variants and the variant_impacts tables.  

============================================================
Using the ``gene_detailed`` & ``gene_summary`` tables
============================================================

---------------------------------------------------------------
Query the gene_detailed table with a join on variants table:
---------------------------------------------------------------
E.g. Get additional transcript info for the most severe impact transcript
e.g. transcript status, transcript start,end and the protein length

.. code-block:: bash

    $ gemini query --header -q "select v.variant_id, v.gene, \
	               v.impact, g.transcript_status, g.transcript, \
				   g.transcript_start, g.transcript_end, g.protein_length, \
				   from variants v, gene_detailed g \
					
				   WHERE v.chrom = g.chrom AND \
						 v.gene = g.gene AND \
						 v.transcript = g.transcript AND \
						 v.impact_severity='HIGH'" test.query.db

	variant_id	gene	impact	transcript_status	transcript	transcript_start	transcript_end	protein_length
	46	SAMD11	frame_shift	KNOWN	ENST00000342066	861118	879955	681
	578	TNFRSF18	frame_shift	PUTATIVE	ENST00000486728	1139224	1141060	169
	733	SCNN1D	stop_gain	NOVEL	ENST00000470022	1217305	1221548	138
	
---------------------------------------------------------------------------
Query the gene_detailed table with a join on the variant_impacts table:
---------------------------------------------------------------------------
E.g. Get the transcript status for all transcripts of the SCNN1D gene where 
impact severity is not 'LOW'.


.. code-block:: bash

    $ gemini query --header -q "select v.gene, g.transcript_status, g.transcript, \
                   v.impact from variant_impacts v, gene_detailed g \
		           
				   WHERE v.transcript = g.transcript AND \
                         v.gene = g.gene AND \
						 v.gene = 'SCNN1D' \
	                     v.impact_severity!='LOW'" test.query.db

	gene	transcript_status	transcript	impact
	SCNN1D	NOVEL	ENST00000470022	non_syn_coding
	SCNN1D	NOVEL	ENST00000470022	frame_shift
	SCNN1D	KNOWN	ENST00000325425	frame_shift
	SCNN1D	KNOWN	ENST00000379116	non_syn_coding
	SCNN1D	KNOWN	ENST00000338555	non_syn_coding
	SCNN1D	KNOWN	ENST00000400928	non_syn_coding
	
---------------------------------------------------------------------------
Query the gene_summary table with a join on the variants table:
---------------------------------------------------------------------------
E.g. Get the synonym/alternate names, RVIS percentile scores and the min-max
start-end of transcripts for genes that have a severely affected transcript 
of a 'HIGH' order.

.. code-block:: bash

    $ gemini query --header -q "select v.chrom, v.gene, g.transcript_min_start, \
	               g.transcript_max_end, g.synonym, g.rvis_pct, v.impact from \
				   variants v, gene_summary g \
                   
				   WHERE v.chrom = g.chrom AND \
                         v.gene = g.gene AND \
                         v.impact_severity='HIGH'" test.query.db

	chrom	gene	transcript_min_start	transcript_max_end	synonym	rvis_pct	impact
	chr1	SAMD11	860260	879955	MGC45873	None	frame_shift
	chr1	TNFRSF18	1138888	1142071	AITR,CD357,GITR	None	frame_shift
	chr1	SCNN1D	1215816	1227409	ENaCdelta,dNaCh	96.77990092	stop_gain
	
-------------------------------------------------------------------------
Query the gene_summary table with a join on the variant_impacts table:
-------------------------------------------------------------------------
E.g. Get all variants of a gene, the affected transcripts and impacts, where 
a mammalian phenotype ID is available for the ``mouse phenotype``.


.. code-block:: bash

    $ gemini query --header -q "select v.variant_id, v.chrom, v.gene, i.impact, \
	               i.transcript, g.mam_phenotype_id from variants v, \
				   variant_impacts i, gene_summary g \
				   
				   WHERE v.variant_id=i.variant_id \
				   i.gene=g.gene AND \
				   v.chrom=g.chrom AND \
				   g.mam_phenotype_id !='None'" test.query.db

	variant_id	chrom	gene	impact	transcript	mam_phenotype_id
	334	chr1	TNFRSF18	non_syn_coding	ENST00000328596	MP:0005397,MP:0005384,MP:0005387
	378	chr1	TNFRSF18	frame_shift	ENST00000486728	MP:0005397,MP:0005384,MP:0005387
	483	chr1	AGRN	synonymous_coding	ENST00000379370	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
	484	chr1	AGRN	exon	ENST00000461111	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
	478	chr1	AGRN	intron	ENST00000461111	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
	479	chr1	AGRN	downstream	ENST00000492947	MP:0005378,MP:0005386,MP:0005388,MP:0005367,MP:0005369,MP:0005371,MP:0003631,MP:0002873,MP:0010768
	
===================================================================
Restrict analysis to transcripts with a valid ``CCDS_ID``
===================================================================
Since the current available transcript sets are more than one (e.g. RefSeq, ENSEMBL and UCSC)
we support information (e.g pathways tool) for the ENSEMBL transcripts but provide a mapping of 
these transcripts to the consensus set agreed upon by all the above three mentioned groups viz. 
``transcripts having a valid CCDS_ID``. Here we show, how we can return variants and their
impacts for only these restricted set of transcripts using the gene_detailed table.


.. code-block:: bash
	
	$ gemini query --header -q "select i.var_id, i.gene, i.impact, i.transcript, \
	                g.transcript_status, ccds_id, g.rvis_pct from \
					variant_impacts i, gene_detailed g where \
					i.transcript=g.transcript and i.gene=g.gene and\
					impact_severity='HIGH' and g.ccds_id!='None'" test.query.db
					
	variant_id	gene	impact	transcript	transcript_status	ccds_id	rvis_pct
	2051	SAMD11	frame_shift	ENST00000342066	KNOWN	CCDS2	None
	3639	CCNL2	splice_acceptor	ENST00000408918	KNOWN	CCDS30558	53.98089172
	3639	CCNL2	splice_acceptor	ENST00000400809	KNOWN	CCDS30557	53.98089172
	13221	SMIM1	frame_shift	ENST00000444870	NOVEL	CCDS57966	None
	21881	NPHP4	splice_acceptor	ENST00000378156	KNOWN	CCDS44052	81.78815758
	
====================================================================
What if I don't see my ``gene`` in the database?
====================================================================
Most genes are popular by their common names while the representation of gene names
in the GEMINI database is mostly HGNC. For e.g ``ARTEMIS`` would be ``DCLRE1C`` in 
the GEMINI database. As such one may miss out on variants if looking for specific
genes by their common names. While, joining the main tables with the gene tables for
synonym information would be useful (as shown in the previous examples), the gene 
tables may also serve as a quick look up for alternate names of a gene, which 
could then be looked up in the database.

.. code-block:: bash
    
	$ gemini query -q "select synonym from gene_summary where \
	                   gene='ARTEMIS'" test.query.db
	A-SCID,FLJ11360,SNM1C,DCLRE1C,SCIDA
	
	or,
	
	$ gemini query -q "select gene from gene_summary where synonym \
	                   like '%ARTEMIS%' and \
	                   is_HGNC='1'" test.query.db
	DCLRE1C
	
	#looking up for DCLRE1C in the database
	$ gemini query -q "select variant_id, chrom, start, end, impact \
	                   from variants where \
	                   gene='DCLRE1C'" test.query.db
 





