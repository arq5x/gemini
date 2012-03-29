---------------------------------------
Stepwise installation and usage of  VEP
---------------------------------------
1. Download the latest version (2.4) of Variant Effect Predictor script from the Ensembl CVS server::

	http://useast.ensembl.org/info/docs/variation/vep/index.html

2. Untar the tarball in the current directory

3. This would create the variant_effect_predictor directory and now do::

	$ cd variant_effect_predictor
	$ perl INSTALL.pl [options]

By default this would install the bioperl-1.2.3, the cache files (in the .vep sub-directory of the users home directory) and the latest version of the Ensembl API (66) (in the variant_effect_predictor directory under a sub-directory 'Bio'). This script is useful for those who do not have all the modules in their system required by VEP, specifically DBI and DBI::mysql. For any alternate options for the installer script, check::

	http://useast.ensembl.org/info/docs/variation/vep/vep_script.html#download
	
4. Users (e.g mac users) who have a problem installing through this script should go for a manual installation of the latest Ensembl API (66) and bioperl-1.2.3 and follow all other installation instructions at::

	http://useast.ensembl.org/info/docs/api/api_installation.html
	
The appropriate pre-build caches should be downloaded to the .vep directory under home from::

	http://useast.ensembl.org/info/docs/variation/vep/vep_script.html#cache
	
To use the cache, the gzip and zcat utilities are required. VEP uses zcat to decompress cached files. For systems where zcat may not be installed or may not work, the following option needs to be added along with the --cache option::

	--compress "gunzip -c"
	
5. You may run the script as::
 
	$ perl variant_effect_predictor.pl [options]
	
We recommend running VEP with the following options as currently we support VEP fields specified as below::

	perl variant_effect_predictor.pl -i example.vcf --cache --compress "gunzip -c" --terms so --condel b --sift b --polyphen b --hgnc --numbers -o output --vcf --fields Consequence,Codons,Amino_acids,Gene,HGNC,Feature,EXON,PolyPhen,SIFT,Condel
	
A documentation of the specified options for VEP may be found at::
	
	http://www.ensembl.org/info/docs/variation/vep/vep_script.html
	
|

------------------------------------------
Stepwise installation and usage of SnpEff
------------------------------------------
To do.
