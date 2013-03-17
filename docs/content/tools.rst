############################
Built-in analysis tools
############################


===========================================================================
``comp_hets``: Identifying potential compound heterozygotes
===========================================================================
Many recessive disorders are caused by compound heterozygotes. Unlike canonical
recessive sites where the same recessive allele is inherited from both parents
at the _same_ site in the gene, compound heterozygotes occur when
the individual's phenotype is caused by two heterozygous recessive alleles at 
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
GTSE1.  The two heterozygotes are reported using the following structure:

``chrom,start,end,ref,alt,genotype,impact,exon,AAF,in_dbsnp``

---------------
``--only_lof``
---------------
By default, all coding variants are explored.  However, one may want to
restrict the analysis to LoF variants using the ``--only_lof`` option.

.. code-block:: bash

	$ gemini comp_hets --only_lof chr22.low.exome.snpeff.100samples.vcf.db

---------------------
``--ignore-phasing``
---------------------
If your genotypes aren't phased, we can't be certain that two heterozygotes
are on opposite alleles.  However, we can still identify pairs of heterozygotes
that are *candidates* for compound heterozygotes. Just use the 
``--ignore-phasing`` option.

----------------------
``--allow-other-hets``
----------------------
By default, the ``comp_hets`` tool will identify candidate pairs of 
heterozygotes that are found in *only one* of the samples in your database. 
Depending on the genetic model, this may be too restrictive.  If you'd like to
identify candidates where other individuals may also be heterozygous, just use
the ``--allow-other-hets`` option

.. code-block:: bash

	$ gemini comp_hets --allow-other-hets chr22.low.exome.snpeff.100samples.vcf.db






===========================================================================
``de_novo``: Identifying potential de novo mutations.
===========================================================================
.. note::

    This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:
    
    ``gemini load -v my.vcf -p my.ped my.db``

To document.





============================================================================
``autosomal_recessive``: Find variants meeting an autosomal recessive model.
============================================================================
.. note::

    This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:
    
    ``gemini load -v my.vcf -p my.ped my.db``

To document.




===========================================================================
``autosomal_dominant``: Find variants meeting an autosomal dominant model.
===========================================================================
.. note::

    This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:
    
    ``gemini load -v my.vcf -p my.ped my.db``

To document.





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
``-t list`` Which regions did a variant overlap?
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
``region``: Extracting variants from specific regions or genes
===========================================================================
One often is concerned with variants found solely in a particular gene or 
genomic region. ``gemini`` allows one to extract variants that fall within 
specific genomic coordinates as follows:

---------
``--reg``
---------
.. code-block:: bash

	$ gemini region --reg chr1:100-200 my.db

----------
``--gene``
----------
Or, one can extract variants based on a specific gene name.

.. code-block:: bash

	$ gemini region --gene PTPN22 my.db


===========================================================================
``windower``: Conducting analyses on genome "windows".
===========================================================================

``gemini`` includes a convenient tool for computing variation metrics across 
genomic windows (both fixed and sliding). Here are a few examples to whet your 
appetite.  If you're still hungry, contact us.

Compute the average nucleotide diversity for all variants found in 
non-overlapping, 50Kb windows.

.. code-block:: bash

	$ gemini windower -w 50000 -s 0 -t nucl_div -o mean my.db

Compute the average nucleotide diversity for all variants found in 50Kb windows 
that overlap by 10kb.

.. code-block:: bash

	$ gemini windower -w 50000 -s 10000 -t nucl_div -o mean my.db




===========================================================================
``db_info``: List the gemini database tables and columns
===========================================================================

Because of the sheer number of annotations that are stored in gemini, there are
admittedly too many columns to remember by rote.  If you can recall the name of 
particular column, just use the ``db_info`` tool.  It will report all of the 
tables and all of the columns / types in each table:

.. code-block:: bash

	$ gemini db_info test.db
	table_name          column_name                   type
	variants            chrom                         text
	variants            start                         integer
	variants            end                           integer
	variants            variant_id                    integer
	variants            anno_id                       integer
	variants            ref                           text
	variants            alt                           text
	variants            qual                          float
	variants            filter                        text
	variants            type                          text
	variants            sub_type                      text
	variants            gts                           blob
	variants            gt_types                      blob
	variants            gt_phases                     blob
	variants            gt_depths                     blob
	variants            call_rate                     float
	variants            in_dbsnp                      bool
	variants            rs_ids                        text
	variants            in_omim                       bool
	variants            clin_sigs                     text
	variants            cyto_band                     text
	variants            rmsk                          text
	variants            in_cpg_island                 bool
	variants            in_segdup                     bool
	variants            is_conserved                  bool
	variants            num_hom_ref                   integer
	variants            num_het                       integer
	variants            num_hom_alt                   integer
	variants            num_unknown                   integer
	variants            aaf                           float
	variants            hwe                           float
	variants            inbreeding_coeff              float
	variants            pi                            float
	variants            recomb_rate                   float
	variants            gene                          text
	variants            transcript                    text
	variants            is_exonic                     bool
	variants            is_coding                     bool
	variants            is_lof                        bool
	variants            exon                          text
	variants            codon_change                  text
	variants            aa_change                     text
	variants            aa_length                     text
	variants            biotype                       text
	variants            impact                        text
	variants            impact_severity               text
	variants            polyphen_pred                 text
	variants            polyphen_score                float
	variants            sift_pred                     text
	variants            sift_score                    float
	variants            anc_allele                    text
	variants            rms_bq                        float
	variants            cigar                         text
	variants            depth                         integer
	variants            strand_bias                   float
	variants            rms_map_qual                  float
	variants            in_hom_run                    integer
	variants            num_mapq_zero                 integer
	variants            num_alleles                   integer
	variants            num_reads_w_dels              float
	variants            haplotype_score               float
	variants            qual_depth                    float
	variants            allele_count                  integer
	variants            allele_bal                    float
	variants            in_hm2                        bool
	variants            in_hm3                        bool
	variants            is_somatic
	variants            in_esp                        bool
	variants            aaf_esp_ea                    float
	variants            aaf_esp_aa                    float
	variants            aaf_esp_all                   float
	variants            exome_chip                    bool
	variants            in_1kg                        bool
	variants            aaf_1kg_amr                   float
	variants            aaf_1kg_asn                   float
	variants            aaf_1kg_afr                   float
	variants            aaf_1kg_eur                   float
	variants            aaf_1kg_all                   float
	variants            grc                           text
	variants            gms_illumina                  float
	variants            gms_solid                     float
	variants            gms_iontorrent                float
	variants            encode_tfbs
	variants            encode_consensus_gm12878      text
	variants            encode_consensus_h1hesc       text
	variants            encode_consensus_helas3       text
	variants            encode_consensus_hepg2        text
	variants            encode_consensus_huvec        text
	variants            encode_consensus_k562         text
	variants            encode_segway_gm12878         text
	variants            encode_segway_h1hesc          text
	variants            encode_segway_helas3          text
	variants            encode_segway_hepg2           text
	variants            encode_segway_huvec           text
	variants            encode_segway_k562            text
	variants            encode_chromhmm_gm12878       text
	variants            encode_chromhmm_h1hesc        text
	variants            encode_chromhmm_helas3        text
	variants            encode_chromhmm_hepg2         text
	variants            encode_chromhmm_huvec         text
	variants            encode_chromhmm_k562          text
	variant_impacts     variant_id                    integer
	variant_impacts     anno_id                       integer
	variant_impacts     gene                          text
	variant_impacts     transcript                    text
	variant_impacts     is_exonic                     bool
	variant_impacts     is_coding                     bool
	variant_impacts     is_lof                        bool
	variant_impacts     exon                          text
	variant_impacts     codon_change                  text
	variant_impacts     aa_change                     text
	variant_impacts     aa_length                     text
	variant_impacts     biotype                       text
	variant_impacts     impact                        text
	variant_impacts     impact_severity               text
	variant_impacts     polyphen_pred                 text
	variant_impacts     polyphen_score                float
	variant_impacts     sift_pred                     text
	variant_impacts     sift_score                    float
	samples             sample_id                     integer
	samples             name                          text
	samples             family_id                     integer
	samples             paternal_id                   integer
	samples             maternal_id                   integer
	samples             sex                           text
	samples             phenotype                     text
	samples             ethnicity                     text
       
