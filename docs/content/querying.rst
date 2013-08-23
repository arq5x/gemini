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

GEMINI has a specific tool for querying a gemini database that has been ``load``ed
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
Filtering on genotypes
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
	 chr1	69510	69511	A	G	OR4F5	./.	./.	A/G	A/G

================================================
Finding out which samples have a variant
================================================
While exploring your data you might hit on a set of interesting variants and want to know
which of your samples have that variant in them. You can display the samples containing
a variant with the --show-sample-variants flag:

.. code-block:: bash

   	$ gemini query --header --show-samples -q "select chrom, start, end, ref, alt \
                                        from variants where is_lof=1 limit 5" test.query.db

	chrom	start	end	ref	alt	variant_samples	HET_samples	HOM_ALT_samples
	chr1	874815	874816	C	CT	1478PC0006B,1478PC0007B,1478PC0010,1478PC0013B,1478PC0022B,1478PC0023B,1478PC0025,1719PC0007,1719PC0009,1719PC0010,1719PC0022	1478PC0006B,1478PC0007B,1478PC0010,1478PC0013B,1478PC0022B,1478PC0023B,1719PC0007,1719PC0009,1719PC0010	1478PC0025,1719PC0022
	chr1	1140811	1140813	TC	T	1478PC0011	1478PC0011
	chr1	1219381	1219382	C	G	1719PC0012	1719PC0012
	chr1	1221487	1221490	CAA	C	1478PC0004	1478PC0004

variant_samples is a list of all of the samples with a variant, HET_samples is the subset
of those heterozygous for the variant and HOM_ALT_samples is the subset homozygous for
the variant.

================================================
--json Reporting query output in JSON format.
================================================
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

To report in JSON format, use the ``--json`` option. For example:

.. code-block:: bash

  $ gemini query --json -q "select chrom, start, end from variants" my.db | head
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



