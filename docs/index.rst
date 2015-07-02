=================================================================
**GEMINI**: *a flexible framework for exploring genome variation*
=================================================================

.. image:: images/overview.png
    :width: 400pt
    :align: center


============================
Latest news
============================
----------------------------
Changes to Inheritance Tools
----------------------------

As of version 0.16.0, the built-in Mendelian inheritance tools are more stringent by default (they can be relaxed with the ``--lenient``) option. By default, samples with unknown phenotype will not affect the tools, and 
strict requirements are placed on family structure. See the the :doc:`docs <content/tools>` for
more info. In addition, the inheritance tools now support multi-generational pedigrees.


----------------------------
Much faster genotype queries
----------------------------
As of version 0.15.0, GEMINI creates auxilliary index files
using the `bcolz` library. This allows queries that filter
on sample genotype information via the ``--gt-filter`` option
to be **up to 1000 times faster than identical queries using
version 0.14.0 or earlier**. Details of the implementation and caveats can be found :doc:`here <content/genotype_query_engines>`.

.. note::

    In order to expedite your queries with the new bcolz indexing strategy, one must invoke the ``--use-bcolz`` 
    option. We will likely phase this option out over time.
    For example::

        gemini query \
          -q "select variant_id, gts.kid from variants"  \
          --gt-filter "gt_types.mom == HET and gt_types.dad == HET and gts.kid != 'HET'" \
          --use-bcolz \
        foo.db


----------------------------
New GEMINI Workflow
----------------------------
As version 0.12.2 of GEMINI it is required that your input VCF file undergo additional preprocessing such that multi-allelic variants are decomposed and normalized using the `vt <http://genome.sph.umich.edu/wiki/Vt>`_ toolset from
the `Abecasis lab <http://genome.sph.umich.edu/wiki/Main_Page>`_. Note that we have also decomposed and normalized all of the VCF-based annotation files (e.g., ExAC, dbSNP, ClinVar, etc.) so that variants and alleles are properly annotated and we minimize false negative and false positive annotations. For a great discussion of why this is necessary, please read `this blog post <http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/>`_ from Eric Minikel in Daniel MacArthur's lab.  

Essentially, VCF preprocessing for GEMINI now boils down to the following steps.

0. If working with GATK VCFs, you need to correct the AD INFO tag definition to play nicely with `vt`.
1. `Decompose <http://genome.sph.umich.edu/wiki/Vt#Decompose>`_ the original VCF such that variants with multiple alleles are expanded into distinct variant records; one record for each REF/ALT combination.
2. `Normalize <http://genome.sph.umich.edu/wiki/Vt#Normalization>`_ the decomposed VCF so that variants are left aligned and represented using the most parsimonious alleles.
3. Annotate with VEP or snpEff.
4. bgzip and tabix.

A workflow for the above steps is given below.

.. code-block:: bash
  
  # setup
  VCF=/path/to/my.vcf
  NORMVCF=/path/to/my.norm.vcf
  REF=/path/to/human.b37.fasta
  SNPEFFJAR=/path/to/snpEff.jar

  # decompose, normalize and annotate VCF with snpEff.
  # NOTE: can also swap snpEff with VEP
  #NOTE: -classic and -formatEff flags needed with snpEff >= v4.1
  zless $VCF \
     | sed 's/ID=AD,Number=./ID=AD,Number=R/'
     | vt decompose -s - \
     | vt normalize -r $REF - \
     | java -Xmx4G -jar $SNPEFFJAR -formatEff -classic GRCh37.75  \
     | bgzip -c > $NORMVCF
  tabix $NORMVCF

  # load the pre-processed VCF into GEMINI
  gemini load --cores 3 -t snpEff -v $NORMVCF $db

  # query away
  gemini query -q "select chrom, start, end, ref, alt, (gts).(*) from variants" \
               --gt-filter "gt_types.mom == HET and \
                            gt_types.dad == HET and \
                            gt_types.kid == HOM_ALT" \
               $db


=================
Overview
=================

GEMINI (GEnome MINIng) is designed to be a flexible framework for exploring genetic variation
in the context of the wealth of genome annotations available for the human genome.
By placing genetic variants, sample genotypes, and useful genome annotations into
an integrated database framework, ``GEMINI`` provides a simple, flexible, yet
very powerful system for exploring genetic variation for disease and
population genetics.

Using the GEMINI framework begins by loading a VCF file (and an optional PED file) into a
database.  Each variant is automatically annotated by comparing it to several
genome annotations from source such as ENCODE tracks, UCSC tracks, OMIM, dbSNP,
KEGG, and HPRD.  All of this information is stored in portable
SQLite database that allows one to explore and interpret both coding and
non-coding variation using "off-the-shelf" tools or an enhanced SQL engine.

Please also see the original `manuscript <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003153>`_.

This `video <http://www.youtube.com/watch?v=p-UWmDG6yj4&t=0m40s>`_ provides more details about GEMINI's aims and
utility.


.. note::

   1. GEMINI solely supports human genetic variation mapped to build 37 (aka hg19) of the human genome.
   2. GEMINI is very strict about adherence to VCF format 4.1.
   3. For best performance, load and query GEMINI databases on the fastest hard drive to which you have access.

=================
Citation
=================
If you use GEMINI in your research, please cite the following manuscript::

    Paila U, Chapman BA, Kirchner R, Quinlan AR (2013) 
    GEMINI: Integrative Exploration of Genetic Variation and Genome Annotations. 
    PLoS Comput Biol 9(7): e1003153. doi:10.1371/journal.pcbi.1003153


=================
Table of contents
=================
.. toctree::
   :maxdepth: 2

   content/installation
   content/quick_start
   content/functional_annotation
   content/preprocessing
   content/querying
   content/tools
   content/browser
   content/database_schema
   content/api
   content/genotype_query_engines
   content/acknowledgements
   content/caveats
   content/history
   content/faq
   content/other
