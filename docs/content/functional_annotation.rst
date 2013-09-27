#############################
Annotation with snpEff or VEP
#############################
.. note::
	
	Version support: VEP versions 2.4 through 73 and core SnpEff
	versions 3.0 through 3.3
	
.. note::
 	
	Version support would be subsequently updated here, as we test along 
	and add or edit changes available with the latest version of these tools.

Stepwise installation and usage of VEP
=======================================

Download the Variant Effect Predictor "standalone perl script" from 
`Ensembl <http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html>`_. 
You can choose a specific version of VEP to download here

Example:

.. code-block:: bash
	
	Download version 71

Untar the tarball into the current directory

.. code-block:: bash
    
    $ tar -zxvf variant_effect_predictor.tar.gz

This will create the variant_effect_predictor directory. Now do the following for install:

.. code-block:: bash

    $ cd variant_effect_predictor
    $ perl INSTALL.pl [options]


By default this would install the API's, bioperl-1.2.3 and the cache files (in 
the $HOME/.vep directory).

Manual installation of VEP
--------------------------

For those (e.g mac users) who have a problem installing through this install script, try a 
manual installation of the API's, BioPerl-1.2.3 and set all pre-requisites for running VEP 
(``DBI`` and ``DBD::mysql`` modules required). The appropriate pre-build caches should be 
downloaded for Human to the ``$HOME/.vep`` directory and then untar.

You may follow instructions at http://www.ensembl.org/info/docs/api/api_installation.html
which provides alternate options for the API installation and additional tips for windows/mac 
users. It also has information for setting up your environment to run VEP.

Example download of the cache files

.. code-block:: bash

	
	$ wget ftp://ftp.ensembl.org/pub/release-73/variation/VEP/homo_sapiens_vep_73.tar.gz
	
You may change the release date in this example to get the appropriate cache files for your version
of VEP that you have installed.

Example

.. code-block:: bash
	
	
	$ wget ftp://ftp.ensembl.org/pub/release-71/variation/VEP/homo_sapiens_vep_71.tar.gz

Cache requires the ``gzip`` and ``zcat`` utilities. VEP uses ``zcat`` 
to decompress cached files. For systems where zcat may not be installed or may 
not work, the following option needs to be added along with the ``--cache`` option:

.. code-block:: bash

    --compress "gunzip -c"

Running VEP
-----------
You may now run VEP as:

.. code-block:: bash

 
    $ perl variant_effect_predictor.pl [OPTIONS]
    
.. note::

We recommend running VEP with the following options as currently we support 
VEP fields specified as below:

.. code-block:: bash

    $ perl variant_effect_predictor.pl -i example.vcf \
       --cache \
       --terms so \
       --sift b \
       --polyphen b \
       --hgnc \
       --numbers \
       -o output \
       --vcf \
       --fields Consequence,Codons,Amino_acids,Gene,HGNC,Feature,EXON,PolyPhen,SIFT
		
	N.B. For VEP version 73, replace "--hgnc" with "--symbol" & "HGNC" in --fields to "SYMBOL"
		
    
A documentation for the above specified options may be found at 
http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html


Stepwise installation and usage of SnpEff
=======================================


.. note::

    Basic Requirements: Java v1.6 or later; at least 2GB of memory


Download the supported versions of SnpEff from http://snpeff.sourceforge.net/download.html 

Example:

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

    $ unzip snpEff_v3_0_GRCh37.66.zip

                                                                                            
To annotate a vcf using snpEff, use the ``default`` options as below:

.. note::

    Memory options for the run may be specified by ``-Xmx2G`` (2GB) or Xmx4G (4GB) 
    based on the requirement

.. code-block:: bash

    $ java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.66 example.vcf > example_snpeff.vcf


If running from a directory different from the installation directory, the 
complete path needs to be specified as,  e.g.:

.. code-block:: bash

    $ java -Xmx4G -jar path/to/snpEff/snpEff.jar -c path/to/snpEff/snpEff.config GRCh37.66 path/to/example.vcf > example_snpeff.vcf


                                                                                                                                                                                                  