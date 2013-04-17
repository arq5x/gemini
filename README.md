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

License
================
`GEMINI` is freely available under the MIT [http://raw.github.com/arq5x/gemini/master/LICENSE](license).