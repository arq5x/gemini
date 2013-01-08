================================================================
**gemini**: *a framework for mining genome variation*
================================================================

=================
Overview
=================
Out overall goal for ``gemini`` is to provide a simple, flexible, and powerful
framework for exploring genetic variation for personalized genomics, as well as
for disease and population genetics.

We accomplish this by integrating genetic variation and genotypes from VCF
files with a wealth of genome annotations from source such as ENCODE,
OMIM, dbSNP, KEGG, and HPRD.  All of this information is stored in portable
SQLite database that allows one to explore and interpret both coding and 
non-coding variation using "off-the-shelf" tools or an enhanced SQL engine.

Caveats and Limitations
=======================
``gemini`` is currently in an *alpha* state.  The basic functionality that we
intend to develop is present, but we may tweak the API and/or restructure 
the underlying database as needed.

Acknowledgements
================
Gemini is developed by Uma Paila and Aaron Quinlan in the 
`Quinlan laboratory <http://quinlanlab.org/>`_ at the University of Virginia.
The development and usability of the project has greatly benefited from input
from Brad Chapman, Oliver Hofmann, and Justin Johnson.

=================
Table of contents
=================
.. toctree::
   :maxdepth: 1

   content/installation
   content/quick_start
   content/functional_annotation
   content/database_schema
   content/bugs
   content/licensing

=============================
Questions and Reporting bugs
=============================
Please contact Aaron Quinlan (aaronquinlan, gmail, .com) 
with questions or issues.
