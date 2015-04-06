GEMINI - integrative exploration of genetic variation and genome annotations.
=============================================================================

Overview
========
The intent of ``GEMINI`` (GEnome MINIing) is to provide a simple, flexible, and 
powerful framework for exploring genetic variation for personal and medical genetics.
GEMINI is unique in that it integrates genetic variation (from VCF files) with
a wealth of genome annotations into a unified database framework. Using this
integrated database as the analysis framework, we aim to leverage the expressive 
power of SQL for data analysis, while attempting to overcome the fundamental 
challenges associated with using databases for very large
(e.g. 1,000,000 variants times 1,000 samples yields one billion genotypes)
datasets. In addition, by defining sample relationships with a PED file, GEMINI allows
one to explore and test for variants that meet specific inheritance models (e.g., 
recessive, dominant, etc.).

Documentation
================
The official documentation is here: http://gemini.readthedocs.org/en/latest/

The following is a video of a high-level talk from SciPy 2013 describing GEMINI.

http://www.youtube.com/watch?v=p-UWmDG6yj4


Citation
================
If you use GEMINI in your research, please cite the following manuscript:

Paila U, Chapman BA, Kirchner R, Quinlan AR (2013) 
GEMINI: Integrative Exploration of Genetic Variation and Genome Annotations. 
PLoS Comput Biol 9(7): e1003153. doi:10.1371/journal.pcbi.1003153

[Manuscript](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003153)


Mailing List
================
Questions and discussion should be directed to the following discussion list:

https://groups.google.com/forum/?fromgroups#!forum/gemini-variation


Acknowledgements
================
GEMINI is being developed in the Quinlan lab (quinlanlab.org) at the University
of Utah and is led by Brent Pedersen, Uma Paila and Aaron Quinlan.  Substantial contributions 
to the code base have also been made by Brad Chapman (@chapmanb) and
Rory Kirchner (@roryk) at the Harvard School of Public Health.


Installation
============
Install ``GEMINI`` using the automated installation script, `gemini_install.py`. This
script installs GEMINI along with required python libraries, third party tools and data 
files used for variant annotation. The installation documentation contains additional 
details on installed files and tools.

http://gemini.readthedocs.org/en/latest/content/installation.html

License
================
`GEMINI` is freely available under the MIT [license](https://github.com/arq5x/gemini/blob/master/LICENSE).


Others
=========
`GEMINI` includes `CADD` scores (PMID: 24487276) for annotating variants.

CADD scores (http://cadd.gs.washington.edu/) are Copyright 2013 University of Washington and 
Hudson-Alpha Institute for Biotechnology (all rights reserved) but are freely available for 
all academic, non-commercial applications. For commercial licensing information contact 
Jennifer McCullar (mccullaj@uw.edu).

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/arq5x/gemini/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

