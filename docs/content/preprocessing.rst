################################################
Preprocessing and Loading a VCF file into GEMINI
################################################

.. _preprocess:

=============================================
Step 1. split, left-align, and trim variants
=============================================

Variants with multiple alternate alleles will not be handled correctly by gemini (or by the tools
used to annotate the variants). As projects get more samples it is likely that a non-negligible 
percentage of site will have multiple alternate alleles.

In addition, variants that are not left-aligned and trimmed can be incorrectly (or not)
annotated.

To reduce the number of false negatives, **we strongly recommend that gemini users split,
left-align, and trim their variants**. The tools we recommend for this are either `vt <https://github.com/atks/vt>`_:

.. code-block:: bash

    vt decompose -s $VCF | vt normalize -r $REFERENCE - > $NEW_VCF

gemini uses the allele depths from the AD tag. In order for `vt` to decompose correctly, users will have
to change the #INFO field for AD in the header from Number=. to Number=R. 

Then the `$NEW_VCF` can be annotated with snpEff or VEP.

===================================
Step 2. Annotate with snpEff or VEP
===================================
.. note::

	Annotate your VCF with SnpEff/VEP, prior to loading it into GEMINI, otherwise the
	gene/transcript features would be set to None.


GEMINI supports gene/transcript level annotations (we do not use pre-computed values here)
from snpEff and VEP  and hence we suggest that you first annotate your VCF with either
of these tools, prior to loading it into GEMINI. The related database columns would be
populated, which would otherwise be set to None if an unannotated VCF file is loaded
into GEMINI.


.. note::
	Choose the annotator as per your requirement!
	Some gene/transcript annotations are available with only one tool (e.g. 
	``Polyphen/Sift`` with VEP). As such these values would be set to None, 
	if an alternate annotator is used during the load step.

Instructions for installing and running these tools can be found in the following section:

:doc:`functional_annotation`

==============================
The basics
==============================

Before we can use GEMINI to explore genetic variation, we must first ``load`` our
VCF file into the GEMINI database framework.  We expect you to have first
annotated the functional consequence of each variant in your VCF using either
VEP or snpEff (Note that v3.0+ of snpEff is required to track the amino acid
length of each impacted transcript). Logically, the loading step is done with
the ``gemini load`` command.  Below are two examples based on a VCF file that
we creatively name my.vcf.  The first example assumes that the VCF has been
pre-annotated with VEP and the second assumes snpEff.

.. code-block:: bash

	# VEP-annotated VCF
	$ gemini load -v my.vcf -t VEP my.db

	# snpEff-annotated VCF
	$ gemini load -v my.vcf -t snpEff my.db

As each variant is loaded into the ``GEMINI`` database framework, it is being
compared against several annotation files that come installed with the software.
We have developed an annotation framework that leverages
`tabix <http://sourceforge.net/projects/samtools/files/tabix/>`_,
`bedtools <http://bedtools.googlecode.com>`_, and
`pybedtools <http://pythonhosted.org/pybedtools/>`_ to make things easy and
fairly performant. The idea is that, by augmenting VCF files with many
informative annotations, and converting the information into a ``sqlite``
database framework, ``GEMINI`` provides a flexible
database-driven API for data exploration, visualization, population genomics
and medical genomics.  We feel that this ability to integrate variation
with the growing wealth of genome annotations is the most compelling aspect of
``GEMINI``.  Combining this with the ability to explore data with SQL
using a database design that can scale to 1000s of individuals (genotypes too!)
makes for a nice, standardized data exploration system.

Many variant callers set filter flags in the VCF file to flag possible problem
variants. By default GEMINI will leave these variants in the database during
loading but they can be filtered out during the loading step by passing the
``--passonly`` flag to load.


================================
Using multiple CPUs for loading
================================

Now, the loading step is very computationally intensive and thus can be very slow
with just a single core.  However, if you have more CPUs in your arsenal,
you can specify more cores.  This provides a roughly linear increase in speed as a
function of the number of cores. On our local machine, we are able to load a
VCF file derived from the exomes of 60 samples in about 10 minutes.  With a
single core, it takes a few hours.


.. note::

    Using multiple cores requires that you have both the ``bgzip`` tool from
    `tabix <http://sourceforge.net/projects/samtools/files/tabix/>`_ and the
    `grabix <https://github.com/arq5x/grabix>`_ tool installed in your PATH.

.. code-block:: bash

    $ gemini load -v my.vcf -t snpEff --cores 20 my.db


=============================================
Using LSF, SGE, SLURM and Torque schedulers
=============================================
One can load VCF files into GEMINI in parallel using many cores on
LSF, SGE, SLURM or Torque clusters. One must simply specify the type of job
scheduler your cluster uses and the queue name to which your jobs
should be submitted.

For example, let's assume you use LSF and a queue named ``preempt_everyone``.
Here is all you need to do:

.. code-block:: bash

    $ gemini load -v my.vcf \
             -t snpEff \
             --cores 50 \
             --queue preempt_everyone \
             --scheduler lsf \
             my.db

===================================
Describing samples with a PED file
===================================
GEMINI also accepts PED files in order to establish the familial relationships
and phenotypic information of the samples in the VCF file.

.. code-block:: bash

    $ gemini load -v my.vcf -p my.ped -t snpEff my.db

The PED file format is documented here: PED_. An example PED file looks like this:

|	1 M10475    None None  1    1
|	1 M10478     M10475  M10500    2    2
|	1 M10500     None    None    2    2
|	1 M128215    M10475  M10500    1    1

The columns are family_id, name, paternal_id, maternal_id, sex and phenotype.

You can also provide a PED file with a heading starting with #, and include extra
fields, like this:

|	#family_id name paternal_id maternal_id sex phenotype hair_color 
| 	1 M10475    None None  1    1 brown 
| 	1 M10478     M10475  M10500    2    2 brown 
| 	1 M10500     None    None    2    2 black 
| 	1 M128215    M10475  M10500    1    1 blue 

This will add the extra columns to the ``samples`` table and allow for you to
use those extra columns during queries.


=======================================
Load GERP base pair conservation scores
=======================================
GERP scores at base pair resolution are loaded by default (Note: This requires a prior install
of the data file by running ``gemini update --dataonly --extra gerp_bp``). However, if not
required, one may optionally skip the load process (to save on the loading time) with the
``--skip-gerp-bp`` option.

.. code-block:: bash

    $ gemini load -v my.vcf --skip-gerp-bp -t snpEff my.db

=========================================
Load CADD scores for deleterious variants
=========================================
CADD scores (http://cadd.gs.washington.edu/) are loaded by default in GEMINI (Note: This requires a
prior install of the data file by running ``gemini update --dataonly --extra cadd_score``). However,
one may optionally skip the load process using the ``--skip-cadd`` option.

.. code-block:: bash

	$ gemini load -v my.vcf --skip-cadd my.db

Updating the samples table in a database
=========================================
If, after loading a database, you find more information about your samples or
want to add a column to the samples table to query on, you can reload the samples
table with a new PED_ file with ``gemini amend --sample``. This is also useful if
you forgot to load a PED_ file when initially loading your database. This file
must have the standard first six columns of a PED_ file, but after that other
columns can be added. The top of the PED_ file also must have a header starting
with # which names all of the columns if there are more than the standard six
PED_ file columns:

.. code-block:: bash

   $ gemini amend --sample your_new_ped_file your.db


===================================
Loading VCFs without genotypes.
===================================
To do.

.. _PED: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
.. _vt-paper: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped 
