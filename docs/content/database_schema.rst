###############
Database Schema
###############


The ``variants`` table
----------------------
================  ========      ====================================================================================
column_name       type          notes
================  ========      ====================================================================================
chrom             STRING        The chromosome on which the variant resides
start             INTEGER       The 0-based start position. 
end               INTEGER       The 1-based end position.
variant_id        INTEGER       PRIMARY_KEY
ref               STRING        Reference allele
alt               STRING        Alternate alele for the variant
qual              INTEGER       Quality score for the assertion made in ALT
filter            STRING        A string of filters passed/failed in variant calling
type              STRING        | The type of variant.
                                | Any of: [*snp*, *indel*]
sub_type          STRING        | The variant sub-type.
                                | If ``type`` is *snp*:   [*ts*, (transition), *tv* (transversion)]
                                | If ``type`` is *indel*: [*ins*, (insertion), *del* (deletion)]
gts               BLOB          A compressed binary vector of sample genotypes (e.g., "A/A", "A|G", "G/G")
gt_types          BLOB          A compressed binary vector of numeric genotype "types" (e.g., 0, 1, 2)
gt_phases         BLOB          A compressed binary vector of sample genotype phases (e.g., False, True, False)
call_rate         FLOAT         The fraction of samples with a valid genotype
in_dbsnp          BOOL          | Is this variant found in dbSnp (build 135)?
                                | 0 : Absence of the variant in dbsnp
                                | 1 : Presence of the variant in dbsnp
rs_ids            STRING        | A comma-separated list of rs ids for variants present in dbsnp
in_omim           BOOL          | 0 : Absence of the variant in OMIM database
                                | 1 : Presence of the variant in OMIM database
clins_sigs        STRING        | A comma-separated list of clinical significance scores for each
                                | of the rs_ids that the variant overlaps in dbsnp. Per dbSNP:
                                | 0 : unknown   1 : untested   2 : non-pathogenic
                                | 3 : probable-non-pathogenic  4 : probable-pathogenic
                                | 5 : pathogenic  6 : drug-response  7 : histocompatibility
                                | 255 : other
cyto_band         STRING        Chromosomal cytobands that a variant overlaps
rmsk              STRING        | A comma-separated list of RepeatMasker annotations that the variant overlaps.
                                | Each hit is of the form: ``name_class_family``
in_cpg_island     BOOL          | Does the variant overlap a CpG island?.
                                | Based on UCSC: Regulation > CpG Islands > cpgIslandExt 
in_segdup         BOOL          | Does the variant overlap a segmental duplication?.
                                | Based on UCSC: Variation&Repeats > Segmental Dups > genomicSuperDups track
num_hom_ref       INTEGER       The total number of of homozygotes for the reference (``ref``) allele
num_het           INTEGER       The total number of heterozygotes observed.
num_hom_alt       INTEGER       The total number of homozygotes for the reference (``alt``) allele
num_unknown       INTEGER       The total number of of unknown genotypes
aaf               FLOAT         The observed allele frequency for the alternate allele
hwe               FLOAT         The Chi-square probability of deviation from HWE (assumes random mating)
pi                FLOAT         The computed nucleotide diversity (pi) for the site
is_exonic         BOOL          Does the variant affect an exon for >= 1transcript?
is_coding         BOOL          Does the variant fall in a coding region (excl. 3' & 5' UTRs) for >= 1 transcript?
is_lof            BOOL          Based on the value of the impact col, is the variant LOF for >= transcript?
depth             INTEGER       The number of aligned sequence reads that led to this variant call
strand_bias       FLOAT         Strand bias at the variant position
rms_map_qual      FLOAT         RMS mapping quality, a measure of variance of quality scores
in_hom_run        INTEGER       Homopolymer runs for the variant allele
num_mapq_zero     INTEGER       Total counts of reads with mapping quality equal to zero
num_alleles       INTEGER       Total number of alleles in called genotypes
num_reads_w_dels  FLOAT         Fraction of reads with spanning deletions
haplotype_score   FLOAT         Consistency of the site with two segregating haplotypes
qual_depth        FLOAT         Variant confidence or quality by depth
allele_count      INTEGER       Allele counts in genotypes
allele_bal        FLOAT         Allele balance for hets
grc               STRING        | Association with patch and fix regions from the Genome Reference Consortium:
                                | http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
                                | Identifies potential problem regions associated with variant calls.
                                | Built with `annotation_provenance/make-ncbi-grc-patches.py`
gms_illumina      FLOAT         | Genome Mappability Scores (GMS) for Illumina error models
                                | Provides low GMS scores (< 25.0 in any technology) from:
                                | http://sourceforge.net/apps/mediawiki/gma-bio/index.php?title=Download_GMS
                                | #Download_GMS_by_Chromosome_and_Sequencing_Technology
                                | Input VCF for annotations prepared with:
                                | https://github.com/chapmanb/bcbio.variation/blob/master/src/bcbio/variation/utils/gms.clj
gms_solid         FLOAT         Genome Mappability Scores with SOLiD error models
gms_iontorrent    FLOAT         Genome Mappability Scores with IonTorrent error models
================  ========      ====================================================================================

|

The ``variant_impacts`` table
----------------------
================  ========      ===============================================================================
column_name       type          notes
================  ========      ===============================================================================
variant_id        INTEGER       PRIMARY_KEY (Foreign key to `variants` table)
anno_id           INTEGER       PRIMARY_KEY (Based on variant transcripts)
gene              STRING        The gene affected by the variant.
transcript        STRING        The transcript affected by the variant.
is_exonic         BOOL          Does the variant affect an exon for this transcript?
is_coding         BOOL          Does the variant fall in a coding region (excludes 3' & 5' UTR's of exons)?
is_lof            BOOL          Based on the value of the impact col, is the variant LOF?
exon              STRING        Exon information for the variants that are exonic
codon_change      STRING        What is the codon change?
aa_change         STRING        What is the amino acid change?
impact            STRING        Impacts due to variation (ref.impact category)
impact_severity   STRING        Severity of the impact based on the impact column value (ref.impact category)
polyphen_pred     STRING        | Impact of the SNP as given by PolyPhen (VEP only) 
                                | benign, possibly_damaging, probably_damaging, unknown
polyphen_scores   FLOAT         Polyphen score reflecting severity (higher the impact, *higher* the score)
sift_pred         STRING        | Impact of the SNP as given by SIFT (VEP only)
                                | neutral, deleterious
sift_scores       FLOAT         SIFT prob. scores reflecting severity (Higher the impact, *lower* the score)
condel_pred       STRING        | Impact of the SNP as given by Condel (VEP only) 
                                | neutral, deleterious
condel_scores     FLOAT         Higher the impact, *higher* the score
================  ========      ===============================================================================

|

The ``samples`` table
----------------------

=============  ==========  ==================================================
column name    type        notes
=============  ==========  ==================================================
sample_id      INTEGER     PRIMARY_KEY
name           STRING      Sample names
family_id      INTEGER     Family ids for the samples [User defined, default: NULL]
paternal_id    INTEGER     Paternal id for the samples [User defined, default: NULL]
maternal_id    INTEGER     Maternal id for the samples [User defined, default: NULL]
sex            STRING      Sex of the sample [User defined, default: NULL]
phenotype      STRING      The associated sample phenotype [User defined, default: NULL]
ethnicity      STRING      The ethnic group to which the sample belongs [User defined, default: NULL]
=============  ==========  ==================================================

|

Details of the ``impact`` and ``impact_severity`` columns
---------------------------------------------------------
================  =======================================
impact severity   impacts
================  =======================================
HIGH              - EXON_DELETED
                  - FRAME_SHIFT
                  - SPLICE_SITE_ACCEPTOR
                  - SPLICE_SITE_DONOR
                  - START_LOST
                  - STOP_GAINED
                  - STOP_LOST
                  - NON_SYNONYMOUS_START
MED               - CODON_CHANGE
                  - CODON_CHANGE_PLUS_CODON_DELETION
                  - CODON_CHANGE_PLUS_CODON_INSERTION
                  - CODON_DELETION
                  - CODON_INSERTION
                  - NON_SYNONYMOUS_CODING
                  - UTR_3_DELETED
                  - UTR_5_DELETED
LOW               - CDS
                  - DOWNSTREAM
                  - EXON
                  - GENE
                  - INTERGENIC
                  - INTERGENIC_CONSERVED
                  - INTRAGENIC
                  - INTRON
                  - INTRON_CONSERVED
                  - START_GAINED
                  - SYNONYMOUS_CODING
                  - SYNONYMOUS_START
                  - SYNONYMOUS_STOP
                  - TRANSCRIPT
                  - UPSTREAM
                  - UTR_3_PRIME
                  - UTR_5_PRIME
================  =======================================



