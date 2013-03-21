
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
Basic Requirements: Java v1.6 or later, memory atleast 2GB


* unix & OS.X: Run the following in a terminal

1. Go to home directory and download the SnpEff version >=3.0

e.g. wget http://sourceforge.net/projects/snpeff/files/snpEff_v3_0_core.zip

N.B. SnpEff should be installed preferably in "snpEff" directory in your home directory. Else, you must update the "data_dir" 
parameter in your snpEff.config file. For e.g. if the installation of snpEff has been done in ~/src instead of ~/ then change 
the data_dir parameter in snpEff.config  as

data_dir = ~/src/snpEff/data/


2. Unzip file
unzip snpEff_v3_0_core.zip


3. Change to snpEff directory and download the genome database
java -jar snpEff.jar download GRCh37.66


4. Unzip the downloaded genome database 
unzip snpEff_v3_2_GRCh37.66.zip

Unzip would create and place the genome in the 'data' directory


5. To annotate a vcf using snpEff do

java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.66 example.vcf > example_snpeff.vcf

N.B. Memory options for the run may be specified by -Xmx2G (2GB) or Xmx4G (4GB) based on the requirement

If running from a directory different from the installation directory, the complete path needs to be specified as:
e.g. java -Xmx4G -jar path/to/snpEff/snpEff.jar -c path/to/snpEff/snpEff.config GRCh37.66 path/to/example.vcf > example_snpeff.vcf



* Windows: Double click on snpEff_v3_0_core.zip and copy contents to the directory of choice. Update the "data_dir" parameter in the snpEff.config file


