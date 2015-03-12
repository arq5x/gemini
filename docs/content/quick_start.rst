###########
Quick start
###########

`gemini` is designed to allow researchers to explore genetic variation contained 
in a `VCF <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_ file.
The basic workflow for working with `gemini` is outlined below.

-------------------------------------------------
Importing VCF files into gemini.
-------------------------------------------------

.. note::

   we now recommend splitting variants with multiple alternates and left-aligning,
   and trimming all variants before loading into gemini.

   See :ref:`preprocess` for a detailed explanation.

Before we can use GEMINI to explore genetic variation, we must first ``load`` our
VCF file into the GEMINI database framework.  We expect you to have first
annotated the functional consequence of each variant in your VCF using either
VEP or snpEff (Note that v3.0+ of snpEff is required to track the amino acid
length of each impacted transcript). Logically,the loading step is done with
the ``gemini load`` command.  Below are two examples based on a VCF file that
we creatively name my.vcf.  The first example assumes that the VCF has been
pre-annotated with VEP and the second assumes snpEff.

.. code-block:: bash

	# VEP-annotated VCF
	$ gemini load -v my.vcf -t VEP my.db

	# snpEff-annotated VCF
	$ gemini load -v my.vcf -t snpEff my.db

Assuming you have a valid VCF file produced by standard variation discovery 
programs (e.g., GATK, FreeBayes, etc.), one loads the VCF into the gemini 
framework with the **load** submodule:

.. code-block:: bash

    $ gemini load -v my.vcf my.db

In this step, `gemini` reads and loads the my.vcf file into a SQLite database 
named my.db, whose structure is described `here <http://gemini.readthedocs.org/en/latest/content/database_schema.html>`_.
While loading the database, `gemini` computes many additional population genetics 
statistics that support downstream analyses. It also stores the genotypes for 
each sample at each variant in an efficient data structure that minimizes the 
database size.

Loading is by far the slowest aspect of GEMINI.  Using multiple CPUs can
greatly speed up this process.

.. code-block:: bash

    $ gemini load -v my.vcf --cores 8 my.db


-------------------------------------------------
Querying the `gemini` database.
-------------------------------------------------

If you are familiar with SQL, ``gemini`` allows you to directly query the database 
in search of interesting variants via the `-q` option.
For example, here is a query to identify all novel, loss-of-function variants 
in your database:

.. code-block:: bash

    $ gemini query -q "select * from variants where is_lof = 1 and in_dbsnp = 0" my.db


Or, we can ask for all variants that substantially deviate from 
Hardy-Weinberg equilibrium:

.. code-block:: bash

    $ gemini query -q "select * from variants where hwe < 0.01" my.db

    
