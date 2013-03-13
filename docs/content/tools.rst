############################
Built-in analysis tools
############################


===========================================================
``annotate``: adding your own custom annotations
===========================================================
It is inevitable that researchers will want to enhance the gemini framework with 
their own, custom annotations. ``gemini`` provides a sub-command called 
``annotate`` for exactly this purpose. As long as you provide a ``tabix``'ed 
annotation file in either BED or VCF format, the annotate tool will, for each 
variant in the variants table, screen for overlaps in your annotation file and 
update a new column in the variants table that you may specify on the command 
line. This is best illustrated by example.

Let's assume you have already created a gemini database of a VCF file using 
the ``load`` module.

.. code-block:: bash

    $ gemini load -v my.vcf -t snpEff my.db
    
Now, let's imagine you have an annotated file in BED format (``crucial.bed``)
that describes regions of the genome that are particularly relevant to your 
lab's research. You would like to annotate in the gemini database which variants 
overlap these crucial regions. We want to store this knowledge in a new column 
in the ``variants`` table called ``crucial_variant`` that tracks whether a given 
variant overlapped (1) or did not overlap (0) intervals in your annotation file.

To do this, you must first TABIX your BED file:

.. code-block:: bash

    $ bgzip crucial.bed
    $ tabix -p bed crucial.bed.gz


------------------------------------------------------
``-t boolean`` Did a variant overlap a region or not?
------------------------------------------------------
Now, you can use this TABIX'ed file to annotate which variants overlap your
crucial regions.  In the example below, the results will be stored in a new
column called "crucial".  The ``-t boolean`` option says that you just want to
track whether (1) or not (0) the variant overlapped one or more of your regions.

.. code-block:: bash

    $ gemini annotate -f crucial.bed.gz -c crucial -t boolean my.db

Since a new columns has been created in the database, we can now directly query 
the new column.  In the example results below, the first and third variants
overlapped a crucial region while the second did not.

.. code-block:: bash

    $ gemini query \
        -q "select chrom, start, end, variant_id, crucial from variants" \
        my.db \
        | head -3
    chr22   100    101    1   1
    chr22   200    201    2   0
    chr22   300    500    3   1
    

-----------------------------------------------------
``-t count`` How many regions did a variant overlap?
-----------------------------------------------------
Instead of a simple yes or no, we can use the ``-t count`` option to *count*
how many crucial regions a variant overlapped.  It turns out that the 3rd
variant actually overlapped two crucial regions.

.. code-block:: bash

    $ gemini annotate -f crucial.bed.gz -c crucial -t count my.db
    
    $ gemini query \
        -q "select chrom, start, end, variant_id, crucial from variants" \
        my.db \
        | head -3
    chr22   100    101    1   1
    chr22   200    201    2   0
    chr22   300    500    3   2

    
-----------------------------------------------------
``-t list`` How many regions did a variant overlap?
-----------------------------------------------------
Lastly, we can *list* which regions a variant overlapped using the ``-t list``
option.  Let's imaging that ``crucial.bed`` looks like this:

.. code-block:: bash
    
    chr22   50    150    crucial1
    chr22   300   400    crucial2
    chr22   350   450    crucial3

When we use ``-t list``, the resulting column can store a comma-separated list
of the region names (column 4).  You can choose whatever column you want to 
store in the database, but in this example, we will use the 4th column (the 
name).  We specify which column to store in the list with the ``-e`` option.

.. code-block:: bash

    $ gemini annotate -f crucial.bed.gz -c crucial -t list -e 4 my.db
    
    $ gemini query \
        -q "select chrom, start, end, variant_id, crucial from variants" \
        my.db \
        | head -3
    chr22   100    101    1   crucial1
    chr22   200    201    2   0
    chr22   300    500    3   crucial2,crucial3


===========================================================================
``comp_hets``: Identifying potential compound heterozygotes
===========================================================================
Many recessive disorders are caused by compound heterozygotes. Unlike canonical
recessive sites where the same recessive allele is inherited from both parents
at the _same_ site in the gene, compound heterozygotes occur when
the individual's phenotype is caused by two heterogeneous recessive alleles at 
_different_ sites in a particular gene. 

So basically, we are looking for two (typically loss-of-function (LoF))
heterozygous variants impacting the same gene at different loci.  The
complicating factor is that this is _recessive_ and as such, we must also 
require that the consequential alleles at each heterozygous site were 
inherited on different chromosomes (one from each parent).  As such, in order
to use this tool, we require that all variants are phased.  Once this has been 
done, the `comp_hets` tool will provide a report of candidate compound
heterozygotes for each sample/gene.

For example:

.. code-block:: bash

	$ gemini comp_hets chr22.low.exome.snpeff.100samples.vcf.db
	sample	gene	het1	het2
	NA19002	GTSE1	chr22,46722400,46722401,G,A,G|A,stop_gain,exon_22,0.005,1	chr22,46704499,46704500,C,A,A|C,stop_gain,exon_22,0.005,0

This indicates that sample NA19002 has a candidate compound heterozygote in
GTSE1.  The two hets are reported using the following structure:

``chrom,start,end,ref,alt,genotype,impact,exon,AAF,in_dbsnp``

By default, all coding variants are explored.  However, one may want to
restrict the analysis to LoF variants.

.. code-block:: bash

	$ gemini comp_hets --only_lof chr22.low.exome.snpeff.100samples.vcf.db


===========================================================================
``region``: Extracting variants from specific regions or genes
===========================================================================
One often is concerned with variants found solely in a particular gene or 
genomic region. ``gemini`` allows one to extract variants that fall within 
specific genomic coordinates as follows:

.. code-block:: bash

	$ gemini region --reg chr1:100-200 my.db

Or, one can extract variants based on a specific gene name.

.. code-block:: bash

	$ gemini region --gene PTPN22 my.db

       