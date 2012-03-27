###########
Quick start
###########

`gemini` is designed to allow researchers to explore genetic variation contained in a `VCF <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_ file.
The basic workflow for working with `gemini` is outlined below.


1. Importing VCF files.
-------------------------------------------------

Assuming you have a valid VCF file produced by standard variation discovery programs (e.g., GATK, FreeBayes, etc.), 
one loads the VCF into the `gemini` framework with the **load** submodule::

    gemini load -v my.vcf my.db

In this step, `gemini` reads and loads the my.vcf file into a SQLite database named my.db, whose structure is described `here <http://nowhere>`_.
While loading the database, `gemini` computes many additional population genetics statistics that support downstream analyses. It also stores
the genotypes for each sample at each variant in an efficient data structure that minimizes the database size.


2. `gemini` queries.
-------------------------------------------------

If you are familiar with SQL, `pop` allows you to directly query the database in search of interesting variants via the `-q` option.
For example, here is a query to identify all novel, loss-of-function variants in your database::

    gemini get -q "select * from variants where is_lof = 1 and in_dbsnp = 0" my.db


Or, we can ask for all variants that substantially deviate from Hardy-Weinberg equilibrium::

    gemini get -q "select * from variants where hwe < 0.01" my.db


3. `gemini` shortcuts.
-------------------------------------------------
We have also developed a number of "shortcuts" (-s option) that allow one to quickly answer common questions without writing 
sophisticated SQL queries::

1. Compute the transition / transversion ratio::

    gemini get -s tstv my.db

2. Principal component analysis::

    gemini get -s pca my.db
    
3. ???::

    gemini get -s ??? my.db

!!!!Uma, add a few more shortcut examples.


4. Custom annotations.
---------------------------------------------------
To do.


    