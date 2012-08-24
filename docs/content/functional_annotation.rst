
Stepwise installation and usage of VEP
---------------------------------------
1. Download the latest version (2.6) of Variant Effect Predictor script from the 
`Ensembl CVS server <http://useast.ensembl.org/info/docs/variation/vep/index.html>`_.

2. Untar the tarball into the current directory.


3. This will create the variant_effect_predictor directory. Now do the following::

    $ cd variant_effect_predictor
    $ perl INSTALL.pl [options]


By default this would install the bioperl-1.2.3, the cache files (in the .vep sub-directory of the users home directory) 
and the latest version of the Ensembl API (68) (in the variant_effect_predictor directory under a sub-directory 'Bio'). 
This script is useful for those who do not have all the modules in their system required by VEP, specifically `DBI` and `DBI::mysql`. 
Use `this <http://useast.ensembl.org/info/docs/variation/vep/vep_script.html#download>`_ link for alternate options of the installer script.


4. Users (e.g mac users) who have a problem installing through this script should go for a manual installation of the latest 
Ensembl API (68) and bioperl-1.2.3 and follow all other installation instructions `here. <http://useast.ensembl.org/info/docs/api/api_installation.html>`_

The appropriate pre-build caches should be downloaded to the `.vep` directory under home from `this <http://useast.ensembl.org/info/docs/variation/vep/vep_script.html#cache>`_ link.

To use the cache, the gzip and zcat utilities are required. VEP uses zcat to decompress cached files. For systems where zcat may not be installed or may not work, the following option needs to be added along with the --cache option::

    --compress "gunzip -c"

5. You may run the script as::
 
    $ perl variant_effect_predictor.pl [OPTIONS]
    
We recommend running VEP with the following options as currently we support VEP fields specified as below::

    perl variant_effect_predictor.pl -i example.vcf \
       --cache --compress "gunzip -c" \
       --terms so \
       --sift b \
       --polyphen b \
       --hgnc \
       --numbers \
       -o output \
       --vcf \
       --fields Consequence,Codons,Amino_acids,Gene,HGNC,Feature,EXON,PolyPhen,SIFT
    
A documentation of the specified options for VEP may be found at http://www.ensembl.org/info/docs/variation/vep/vep_script.html


Stepwise installation and usage of SnpEff
------------------------------------------
To do.
