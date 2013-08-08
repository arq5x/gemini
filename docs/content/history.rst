#############################
Release History
#############################


0.5.0b (2013-Jul-23)
=======================================
1. Tolerate either -9 or 0 for unknown parent or affected status in PED files.
2. Refine the rules for inheritance and parental affected status for autosomal dominant inheritance models.
3. The ``autosomal_dominant``, ``autosomal_recessive``, and ``de_novo`` mutation tools have received the following improvements.

    -  improved speed (especially when there are multiple families)
    -  by default, all columns in the variant table are reported and no conditions
       are placed on the returned variants.  That is, as long as the variant meets 
       the inheritance model, it will be reported.
    -  the addition of a ``--columns`` option allowing one to override the above default
       behavior and report a subset of columns.
    -  the addition of a ``--filter`` option allowing one to override the above default
       behavior and filter reported variants based on specific criteria.

4. The default minimum aligned sequencing depth for each variant reported by 
the ``de_novo`` tool is 0.  Greater stringency can be applied with the ``-d``
option.

0.4.0b (2013-Jun-12)
=======================================
1. Added new ``gt_ref_depths``, ``gt_alt_depths``, and ``gt_quals`` columns.
2. Added a new ``--show-samples`` option to the ``query`` module to display samples with alternate allele genotypes.
3. Improvements and bug fixes for installation.

0.3.0b
=======================================
1. Improved speed for adding custom annotations.
2. Added GERP conserved elements.
3. Optionally addition of GERP conservation scores at base pair resolution.
4. Move annotation files to Amazon S3.

