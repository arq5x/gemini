GEMINI - integrative exploration of genetic variation and genome annotations.
=============================================================================

Overview
========
The intent of ``GEMINI`` (GEnome MINIing) is to provide a simple, flexible, and 
powerful framework for exploring genetic variation for personal and medical genetics.
GEMINI is unique in that it integrated genetic variation (from VCF files) with
a wealth of genome annotations into a unified database framework. Using this
integrated database as the analysis framework, we aim to leverage the expressive 
power of SQL for data analysis, while attempting to overcome the fundamental 
challenges associated with using databases for very large
(e.g. 1,000,000 variants times 1,000 samples yields one billion genotypes)
datasets.

Documentation
================
The official documentation is here: http://gemini.readthedocs.org/en/latest/

Mailing List:
https://groups.google.com/forum/?fromgroups#!forum/gemini-variation

Acknowledgements
================
GEMINI is being developed in the Quinlan lab (quinlanlab.org) at the University
of Virginia and is led by Uma Paila and Aaron Quinlan.  Substantial contributions 
to the code base have also been made by Brad Chapman (@chapmanb) and
Rory Kirchner (@roryk) at the Harvard School of Public Health.


Installation
============
Install ``GEMINI`` using the automated installation script, `gemini_install.py`. This
script installs GEMINI along with required python libraries, third party tools and data 
files used for variant annotation. The installation documentation contains additional 
details on installed files and tools.::

    wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
    python gemini_install.py /usr/local /usr/local/share/gemini

If you don't have administrative (or sudo) priveleges on you machine, use the
`--nosudo` option::

    python gemini_install.py --nosudo ~/src ~/data/gemini

If using the `--nosudo` option, you will need to add the `gemini/bin` directory
that was created as part of the installation to your PATH.

The first directory specified (`/usr/local` in the first example above) dictates 
the directory in which the gemini software should be installed.  The second
directory specified (`/usr/local/share/gemini` in the first example above) specifies
the directory in which the gemini annotations files should be installed.  NOTE that
these are rather large files (some > 1Gb).

The [installation documentation][install] has additional tips and tricks if you
run into any problems downloading or using the installation script.

[install]: http://gemini.readthedocs.org/en/latest/content/installation.html

License
================
`GEMINI` is freely available under the MIT [license](https://github.com/arq5x/gemini/blob/master/LICENSE).
