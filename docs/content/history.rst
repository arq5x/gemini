#############################
Release History
#############################

(Future)
=======================================
#. Support for both SO and HGVS terms.
#. Update parsing for snpEff >= 4.0
#. Speed up and generalize database loading for multiple genome builds and species.
#. Add an `is_splicing` column.

0.16.0
=======================================
1. The built-in inheritance model tools (``auto_rec``, etc.) have been modified to be more
   restrictive in order to remove false positive candidates. The strictness can be reduced by using the ``--lenient`` option.
2. Leverage bcolz indexing for the built-in inheritance model tools to dramatically improve speed.
3. Support for multi-generational pedigrees for the built in inheritance model tools. (thanks to Jessica, Andrew,
   and jmcelwee for extensive discussion https://github.com/arq5x/gemini/issues/388)
4. Leverage genotype likelihoods in tools other than ``mendel_errors`` as a means to filter variants.
5. Automatically phase genotypes by transmission on the fly for the `comp_hets` tool.
6. Further performance improvements for bcolz queries
7. The ``--affected-only`` option has been made the default and it's opposing replacement named ``--allow-unaffected`` to revert.
8. Fixed a reporting error for the inheritance tools (i.e., family_id was mis-specified in output).
9. Annotate the variants table with impact even if there is not severe impact. Thanks to @mjsduncan for reporting.
10. Reduce memory requirements when loading. Thanks to @mjsduncan for reporting.

0.15.1
======
1. Fix regression in grabix. Thanks to Sven-Eric Shelhorn for reporting.
2. Fix handling of samples with "-". Thanks to Uma Paila for reporting.

0.15.0
=======================================
1. Use external index to speed genotype queries (this is created by default on load unless --no-bcolz is specified)
2. Match on ref and alternate alleles (not just position) when annotating with VCF. Thanks Jeremy Goecks.
3. Related to matching, we now load extra annotation, e.g. VEP as VCF and require ref and alt matching. Previously was done with bed overlap.
4. Faster queries due to lazy loading of genotype columns.
5. Read gt\* columns from the database for better backward compatibility.
6. Code cleanup. Thanks to Christian Brueffer.

0.14.0 
=======================================
1. Standardized the output from the built-in tools into a common, BED+ format. Thanks to feedback from Jessica Chong and Daniel Gaston.
2. Release of `mendel_errors` tool which also outputs the type of error and the probability (based on PL's)
3. Improvements to the `load` tool when running on large compute clusters using PBS, SGE, SLURM, etc. Also provde a workaround for NFS locking issues. Many thanks to Ben Weisburd in Daniel Macarthur's lab.
4. Improve preprocess script to support varscan, platypus (https://gist.github.com/brentp/4db670df147cbd5a2b32)
5. Performance improvements for many of the built-in tools (pre-compile evals)
6. Bug fix for installation with sudo privileges.


0.13.1 (2015-Apr-09)
=======================================
1. Major `query` speed improvements. For example, the following query goes from 43 seconds in version 0.12.2 to 11 seconds in 0.13.0. All queries involving `gt_*` fields should be substantially faster.
  ::

    $ gemini query \
            -q "select chrom, start, (gts).(*) from variants" data/tmaster.db \
            --gt-filter "(gt_depths).(*).(>=20).(all)" > /dev/null

2. Speed improvements to `load`. The following went from 7 minutes 9 seconds to 6 minutes 21 seconds.
  ::

    $ gemini load -t VEP -v data/v100K.vcf.gz data/tmaster.db --cores 4

3. We added the `gt_phred_ll_homref`, `gt_phred_ll_het`, `gt_phred_ll_homalt` columns to database. These are the genotype likelihoods pulled from the GL or PL columns of the VCF if available. They can all be queried and filtered in the same way as existing gt_* columns. In future releases, we are planning tp use genotype likelihood to assign likelihoods to de novo mutations, mendelian violations, and variants meeting other inheritance patterns.

4. Fixed bugs related to splitting multiple alts (thanks to @jdh237)

5. We are working to improve development and release testing. This is ongoing, but we now support gemini_install.py --version unstable so that users can try out the latest changes and help with testing before releases. gemini_update is still limited to master as the most recent version.

6. Update cyvcf so it doesn't error when AD tag is used for non-list data.

7. Fix regression in cyvcf to handle Flags in info field. (Thanks to Jon for reporting)

8. Improvements to install related to PYTHONHOME and other env variables(@chapmanb & @bw2)



0.12.2 
=======================================
Corrected a stale .c file in the cyvcf library. This is effectively a replacement for the 0.12.1 release.


0.12.1 
=======================================
1. Support for input VCF files containing variants with multiple alternate alleles. Thanks to Brent Pedersen.
2. Updated, decomposed, and normalized the ExAC, Clinvar, Cosmic, dbSNP, and ESP annotation files to properly support variants with multiple alternate alleles.
3. Updated the logic for the clinvar significance column to retain all documented singificances.
4. Support for VCF annotation files in the `annotate` tool.
5. Improved the speed of loading by 10-15%. Thanks to Brent Pedersen.
6. Added `--only-affected` and `--min-kindreds` options to the compound heterozygotes tool.
7. Added a `--format vcf` option to the `query` tool to output query results in VCF format.
8. Added the `--families` option to the `auto_*`, `de_novo`, and `comp_hets` tools. Thanks to Mark Cowley and Tony Roscioli.
9. Added the `--only-affected` option to the `de_novo` tool.
10. Allow the `--sample-filter` to work with `--format TPED`. Thanks to Rory Kirchner.
11. Add `--format sampledetail` option that provides a melted/tidy/flattened version of samples along with `--showsample` and includes information from samples table. Thanks to Brad Chapman.
12. Add 'not' option to --in filtering. Thanks to Rory Kirchner.
13. Fixed a bug in the `de_novo` tool that prevented proper function when families have affected and unaffected children. Thanks to Andrew Oler.
14. Fixed a bug in cyvcf that falsely treated '.|.' genotypes as homozygous alternate.  Thanks to Xiao Xu.
15. GEMINI now checks for and warns of old grabix index files. Thanks to Andrew Oler and Brent Pedersen.
16. Fixed a bug that added newlines at the end of tab delimited PED files.  Thanks to Brad Chapman.


0.11.0
=======================================
1. Integration of ExAC annotations (v0.2): http://exac.broadinstitute.org/
2. New tools for cancer genome analysis. Many thanks to fantastic work from Colby Chiang.
  - `gemini set_somatic`
  - `gemini actionable_mutations`
  - `gemini fusions`
3. Improved support for structural variants. New columns include:
  - `sv_cipos_start_left`
  - `sv_cipos_end_left`
  - `sv_cipos_start_right`
  - `sv_cipos_end_right`
  - `sv_length`
  - `sv_is_precise`
  - `sv_tool`
  - `sv_evidence_type`
  - `sv_event_id`
  - `sv_mate_id`
  - `sv_strand`
4. Updated the 1000 Genomes annotations to the Phase variant set.
5. Added `clinvar_causal_allele` column.
6. Fixed a bug in grabix that caused occasional duplicate and missed variants.

0.10.1
=======================================
1. Add `fitCons <http://biorxiv.org/content/early/2014/09/11/006825>` scores as
   an additional measure of potential function in variants of interest, supplementing
   existing CADD and dbNSFP approaches.
2. Updated Clinvar, COSMIC, and dbSNP to their latest versions.


0.10.0 
===================
1. Provide an ``--annotation-dir`` argument that specifies the path the
   annotation databases, to overwrite configured data inputs. Thanks to Björn Grüning,
2. Support reproducible versioned installs of GEMINI with Python
   dependencies. Enables Galaxy integration. Thanks to Björn Grüning,


0.8.0
=======================================
1. Support arbitrary annotation supplied to VEP, which translate into queryable
   columns in the main variant table.
2. Improve the power of the genotype filter wildcard functionality.


0.7.1
=======================================
1. Suppress openpyxl/pandas warnings (thanks to @chapmanb)
2. Fix unit tests to account for cases where a user has not downloaded the CADD or GERP annotation files.
   Thanks to Xialoin Zhu and Daniel Swensson for reporting this and to Uma Paila for correcting it.

0.7.0
=======================================
1. Added support for CADD scores via new ``cadd_raw`` and ``cadd_scaled`` columns.
2. Added support for genotype wildcards in query select statements. E.g., ``SELECT chrom, start, end (gts).(phenotype==2) FROM variants``. See http://gemini.readthedocs.org/en/latest/content/querying.html#selecting-sample-genotypes-based-on-wildcards.
3. Added support for genotype wildcards in the --gt-filter. E.g., ``--gt-filter "(gt_types).(phenotype==2).(==HET)``. See http://gemini.readthedocs.org/en/latest/content/querying.html#gt-filter-wildcard-filtering-on-genotype-columns.
4. Added support for the VCF INFO field both in the API and as a column that can be SELECT'ed.
5. Upgraded to the latest version of ClinVar.
6. Standardized impacts to use Sequence Ontology (SO) terms.
7. Automatically add indexes to custom, user-supplied annotation columns.
8. Improvements to the installation script.
9. Fixed bugs in the handling of ClinVar UTF8 encoded strings.
10. Upgraded the ``gene_summary`` and ``gene_detailed`` tables to version 75 of Ensembl.
11. Added support for the MPI Mouse Phenotype database via the ``mam_phenotype_id`` column in the ``gene_summary`` table.
12. Enhanced security.
13. Corrected the ESP allele frequencies to be based report _alternate_ allele frequency instead of _minor_ allele frequency.
14. VEP version support updated (73-75) Support for aa length and bio type in VEP. 
15. The `lof_sieve` tool support has been extended to VEP annotations.
16. Added the ``ccds_id`` and ``entrez_id`` columns to the ``gene_detailed`` table.


0.6.6
=======================================
1. Added COSMIC mutation information via new cosmic_ids column.


0.6.4 (2014-Jan-03)
=======================================

1. New annotation: experimentally validated human enhancers from VISTA.
2. Installation improvements to enable isolated installations inside of virtual
   machines and containers without data. Allow data-only upgrades as part of
   ``update`` process.
3. Fix for gemini query error when ``--header`` specified (#241).

0.6.3.2 (2013-Dec-10)
=======================================
1. Fixed a bug that caused ``--gt-filter`` to no be enforced from ``query`` tool unless a GT* column was selected.
2. Support for ref and alt allele depths provided by FreeBayes.

0.6.3.1 (2013-Nov-19)
=======================================
1. Fixed undetected bug preventing the ``comp_hets`` tool from functioning.
2. Added unit tests for the ``comp_hets`` tool.

0.6.3 (2013-Nov-7)
=======================================
1. Addition permutation testing to the c-alpha test via the ``--permutations``
   option.
2. Addition of the ``--passonly`` option during loading to filter out all
   variants with a filter flag set.
3. Fixed bug with parallel loading using the extended sample table format.
4. SLURM support added.
5. Refactor of loading options to remove explosion of xxx-queue options. Now
   load using ``--scheduler`` on ``--queue``.
6. Refactor of Sample class to handle the expanded samples table.
7. Addition of ``--carrier-summary-by-phenotype`` for summarizing the counts of
   carriers and non-carriers stratified by the given sample phenotype column.
8. Added a ``--nonsynonymous`` option to the C-alpha test.
9. Added ``gemini amend`` to edit an existing database. For now only handles updating
   the samples table.
10. Fixed a bug that prevented variants that overlapped with multiple 1000G variants
    from having AAF info extracted from 1000G annotations.  This is now corrected such
    that multiple overlaps with 1000G variants are tolerated, yet the logic ensures 
    that the AAF info is extracted for the correct variant.
11. Fixed installation issues for the GEMINI browser.
12. Added ``--show-families`` option to gemini query.


0.6.2 (2013-Oct-7)
=======================================
1. Moved `--tped` and `--json` options into the more generic `--format` option.
2. Fixed bug in handling missing phenotypes in the sample table.
3. Fixed `--tped` output formatting error.
4. API change: GeminiQuery.run takes an optional list of predicates that a row
   must pass to be returned.
5. `--sample-filter` option added to allow for restricting variants to samples
   that pass the given sample query.
6. ethnicity removed as a default PED field.
7. PED file format extended to allow for extra columns to be added to the samples table under the column named in the header.
8. The autosomal_recessive and autosomal_dominant tools now warn, but allow for variants to be detected in the absence of known parent/child relationships.


0.6.1 (2013-Sep-09)
=======================================
1. Corrected bug in de_novo tool that was undetected in 0.6.0.  Unit tests have been added to head this off in the future. Thanks to **Jessica Chong**
2. Added the `-d` option (minimum sequence depth allowed for a genotype) to the `autosmal_recessive` and `autosmal_dominant` tools.
3. New `--tped` option in the `query` tool for reporting variants in TPED format. Thanks to **Rory Kirchner**.
4. New `--tfam` option in the `dump` tool for reporting sample infor in TFAM format. Thanks to **Rory Kirchner**.



0.6.0 (2013-Sep-02)
=======================================
1. Add the ``--min-kindreds`` option to the ``autosomal_recessive`` and ``autosomal_dominant`` tools to restrict candidate variants/genes to those affecting at least ``--min-kindreds``. Thanks to **Jessica Chong**
2. Addition of a new ``burden`` tool for gene or region based burden tests.  First release supports the C-alpha test.  Thanks to **Rory Kirchner**.
3. Use of Continuum Analytics Anaconda python package for the automated installer. Thanks to **Brad Chapman**.
4. Enhancements to the ``annotate`` tool allowing one to create new database columns from values in custom BED+ annotation files.  Thanks to **Jessica Chong** and **Graham Ritchie**.
5. Addition of the ``--column``, ``--filter``, and ``--json`` options to the ``region`` tool.
6. Improvements to unit tests.
7. Allow alternate sample delimiters in the ``query`` tool via the ``--sample-delim`` option.  Thanks to **Jessica Chong**.
8. Provide a REST-like interface to the gemini browser.  In support of future visualization tools.
9. Allow the ``query`` tool to report results in JSON format via the ``--json`` option.
10. Various minor improvements and bug fixes.




0.5.0b (2013-Jul-23)
=======================================
1. Tolerate either -9 or 0 for unknown parent or affected status in PED files.
2. Refine the rules for inheritance and parental affected status for autosomal dominant inheritance models.
3. The ``autosomal_dominant``, ``autosomal_recessive``, and ``de_novo`` mutation tools have received the following improvements.

    -  improved speed (especially when there are multiple families)
    -  by default, all columns in the variant table are reported and no conditions are placed on the returned variants.  That is, as long as the variant meets
       the inheritance model, it will be reported.
    -  the addition of a ``--columns`` option allowing one to override the above default behavior and report a subset of columns.
    -  the addition of a ``--filter`` option allowing one to override the above default behavior and filter reported variants based on specific criteria.

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

