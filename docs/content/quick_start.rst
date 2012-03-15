###########
Quick start
###########

`pop` is designed to allow researchers to explore genetic variation contained in a `VCF <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_ file.
The basic workflow for working with pop is outlined below.


1. Importing a VCF file into the `pop` framework.
-------------------------------------------------

Assuming you have a valid VCF file produced by standard variation discovery programs (e.g., GATK, FreeBayes, etc.), 
one loads the VCF into the `pop` framework with the **load** submodule::

    pop load -v my.vcf my.db

In this step, `pop` reads and loads the my.vcf file into a SQLite database named my.db, whose structure is described `here <http://nowhere>`_.
While loading the database, `pop` computes many additional population genetics statistics that support downstream analyses. It also stores
the genotypes for each sample at each variant in an efficient data structure that minimizes the database size.


2. Querying the `pop` database using SQL queries.
-------------------------------------------------

If you are familiar with SQL, `pop` allows you to directly query the database in search of interesting variants via the `-q` option.
For example, here is a query to identify all novel, loss-of-function variants in your database::

    pop get -q "select * from variants where is_lof = 1 and in_dbsnp = 0" my.db


Or, we can ask for all variants that substantially deviate from Hardy-Weinberg equilibrium::

    pop get -q "select * from variants where hwe < 0.01" my.db


3. Querying the `pop` database using "shortcuts".
-------------------------------------------------
We have also developed a number of "shortcuts" (-s option) that allow one to quickly answer common questions without writing 
sophisticated SQL queries::

1. Compute the transition / transversion ratio::

    pop get -s tstv my.db

2. Principal component analysis::

    pop get -s pca my.db
    
3. ???::

    pop get -s ??? my.db

!!!!Uma, add a few more shortcut examples.


4. Adding your own annotations to the database.
---------------------------------------------------
To do.


    