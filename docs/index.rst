=================================================================
**gemini**: *a flexible framework for exploring genome variation*
=================================================================

=================
Overview
=================

Gemini is designed to be a flexible framework for exploring genetic variation
in the context of the wealth of genome annotations available for the human genome.
By placing genetic variants, sample genotypes, and useful genome annotations into
an integrated database framework, ``gemini`` provides a simple, flexible, yet
very powerful system for exploring genetic variation for for disease and 
population genetics.

Using the Gemini framework begins by loading a VCF file into a
database.  Each variant is automatically annotated by comparing it to several
genome annotations from source such as ENCODE tracks, UCSC tracks, OMIM, dbSNP, 
KEGG, and HPRD.  All of this information is stored in portable
SQLite database that allows one to explore and interpret both coding and 
non-coding variation using "off-the-shelf" tools or an enhanced SQL engine.

=================
Table of contents
=================
.. toctree::
   :maxdepth: 2

   content/installation
   content/quick_start
   content/loading
   content/querying
   content/tools
   content/browser
   content/database_schema
   content/api
   content/acknowledgements
   content/caveats
   content/bugs
   content/licensing


.. content/functional_annotation
