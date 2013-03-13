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


       