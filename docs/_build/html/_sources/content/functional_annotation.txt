#############################
Annotation with snpEff or VEP
#############################


Stepwise installation and usage of VEP
=======================================

Download the latest version of Variant Effect Predictor "standalone Perl script"
from the `Ensembl CVS server <http://useast.ensembl.org/info/docs/variation/vep/index.html>`_.
For example:

.. code-block:: bash

    $ open http://useast.ensembl.org/info/docs/variation/vep/index.html
    
Untar the tarball into the current directory.

.. code-block:: bash
    
    $ tar -zxvf variant_effect_predictor.tar.gz

This will create the variant_effect_predictor directory. Now do the following:

.. code-block:: bash

    $ cd variant_effect_predictor
    $ perl INSTALL.pl [options]


By default this would install the bioperl-1.2.3, the cache files (in the .vep 
sub-directory of the users home directory) and the latest version of the 
Ensembl API (68) (in the ``variant_effect_predictor`` directory under a 
sub-directory ``Bio``). This script is useful for those who do not have all the 
modules in their system required by VEP, specifically `DBI` and `DBI::mysql`. 
Use `this <http://useast.ensembl.org/info/docs/variation/vep/vep_script.html#download>`_ 
link for alternate options of the installer script.


Users (e.g mac users) who have a problem installing through this script 
should go for a manual installation of the latest Ensembl API (68) and 
bioperl-1.2.3 and follow all other installation instructions 
`here. <http://useast.ensembl.org/info/docs/api/api_installation.html>`_

The appropriate pre-build caches should be downloaded to the `.vep` directory 
under home from `this <http://useast.ensembl.org/info/docs/variation/vep/vep_script.html#cache>`_ link.

To use the cache, the gzip and ``zcat`` utilities are required. VEP uses ``zcat`` 
to decompress cached files. For systems where zcat may not be installed or may 
not work, the following option needs to be added along with the ``--cache`` option:

.. code-block:: bash

    --compress "gunzip -c"

You may run the script as:

.. code-block:: bash

 
    $ perl variant_effect_predictor.pl [OPTIONS]
    
.. note::

We recommend running VEP with the following options as currently we support 
VEP fields specified as below:

.. code-block:: bash

    $ perl variant_effect_predictor.pl -i example.vcf \
       --cache --compress "gunzip -c" \
       --terms so \
       --sift b \
       --polyphen b \
       --hgnc \
       --numbers \
       -o output \
       --vcf \
       --fields Consequence,Codons,Amino_acids,Gene,HGNC,Feature,EXON,PolyPhen,SIFT
    
A documentation of the specified options for VEP may be found at 
http://www.ensembl.org/info/docs/variation/vep/vep_script.html


Stepwise installation and usage of SnpEff
=======================================


.. note::

    Basic Requirements: Java v1.6 or later; at least 2GB of memory


Go to home directory and download the SnpEff version >=3.0. For example:

.. code-block:: bash

    $ wget http://sourceforge.net/projects/snpeff/files/snpEff_v3_0_core.zip

.. note::

    SnpEff should be installed preferably in ``snpEff`` directory in your 
    home directory. Else, you must update the ``data_dir`` parameter in 
    your snpEff.config file. For e.g. if the installation of snpEff has been done 
    in ``~/src`` instead of ``~/`` then change the data_dir parameter in 
    snpEff.config to ``data_dir = ~/src/snpEff/data/``


Unzip the downloaded package.

.. code-block:: bash

    $ unzip snpEff_v3_0_core.zip


Change to the ``snpEff`` directory and download the genome database.

.. code-block:: bash

    $ cd snpEff_v3_0_core
    $ java -jar snpEff.jar download GRCh37.66


Unzip the downloaded genome database. This will create and place the genome 
in the 'data' directory

.. code-block:: bash

    $ unzip snpEff_v3_2_GRCh37.66.zip


To annotate a vcf using snpEff, use the following command:

.. note::

    Memory options for the run may be specified by ``-Xmx2G`` (2GB) or Xmx4G (4GB) 
    based on the requirement

.. code-block:: bash

    $ java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.66 example.vcf > example_snpeff.vcf


If running from a directory different from the installation directory, the 
complete path needs to be specified as,  e.g.:

.. code-block:: bash

    $ java -Xmx4G -jar path/to/snpEff/snpEff.jar -c path/to/snpEff/snpEff.config GRCh37.66 path/to/example.vcf > example_snpeff.vcf


