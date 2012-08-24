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
anno_id           INTEGER       Variant transcript number for the most severely affected transcript
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
is_conserved      BOOL          | Does the variant overlap a conserved region?
                                | Based on the 29-way mammalian conservation study
num_hom_ref       INTEGER       The total number of of homozygotes for the reference (``ref``) allele
num_het           INTEGER       The total number of heterozygotes observed.
num_hom_alt       INTEGER       The total number of homozygotes for the reference (``alt``) allele
num_unknown       INTEGER       The total number of of unknown genotypes
aaf               FLOAT         The observed allele frequency for the alternate allele
hwe               FLOAT         The Chi-square probability of deviation from HWE (assumes random mating)
inbreeding_coeff  FLOAT         The inbreeding co-efficient that expresses the likelihood of effects due to inbreeding
pi                FLOAT         The computed nucleotide diversity (pi) for the site
recomb_rate       FLOAT         | Returns the mean recombination rate at the variant site
                                | Based on HapMapII_GRCh37 genetic map
gene              STRING        Corresponding gene name of the highly affected transcript
transcript        STRING        The variant transcript that was most severely affected
is_exonic         BOOL          Does the variant affect an exon for >= 1transcript?
is_coding         BOOL          Does the variant fall in a coding region (excl. 3' & 5' UTRs) for >= 1 transcript?
is_lof            BOOL          Based on the value of the impact col, is the variant LOF for >= transcript?
exon              STRING        Exon information for the severely affected transcript
codon_change      STRING        What is the codon change?
aa_change         STRING        What is the amino acid change (for an snp)?
aa_length         STRING        The length of CDS in terms of number of amino acids
biotype           STRING        The 'type' of the severely affected transcript (e.g.protein-coding, pseudogene, rRNA etc.)
impact            STRING        The consequence of the most severely affected transcript
impact_severity   STRING        Severity of the highest order observed for the variant
polyphen_pred     STRING        Polyphen predictions for the snps (only with VEP) for the severely affected transcript
polyphen_score    FLOAT         Polyphen scores for the severely affected transcript
sift_pred         STRING        SIFT predictions for the snp's (VEP only) for the most severely affected transcript
sift_score        FLOAT         SIFT scores for the predictions
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
in_esp            BOOL          Presence/absence of the variant in the ESP project data
aaf_esp_ea        FLOAT         Minor Allele Frequency of the variant for European Americans in the ESP project
aaf_esp_aa        FLOAT         Minor Allele Frequency of the variant for African Americans in the ESP project
aaf_esp_all       FLOAT         Minor Allele Frequency of the variant w.r.t both groups in the ESP project
exome_chip        BOOL          Whether an SNP is on the Illumina HumanExome Chip
in_1kg            BOOL          Presence/absence of the variant in the 1000 genome project data
aaf_1kg_amr       FLOAT         Allele Frequency of the variant for samples in AMR based on AC/AN (1000g project)
aaf_1kg_asn       FLOAT         Allele frequency of the variant for samples in ASN based on AC/AN (1000g project)
aaf_1kg_afr       FLOAT         Allele frequency of the variant for samples in AFR based on AC/AN (1000g project)
aaf_1kg_eur       FLOAT         Allele Frequency of the variant for samples in EUR based on AC/AN (1000g project)
aaf_1kg_all       FLOAT         Global allele frequency (based on AC/AN) (1000g project) 
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
encode_tfbs       STRING        | Comma-separated list of transcription factors that were
                                | observed by ENCODE to bind DNA in this region.  Each hit in the list is constructed
                                | as TF_MAXSCORE_CELLCOUNT, where:
                                |   *TF* is the transcription factor name
                                |   *MAXSCORE* is the highest signal strength observed in any of the cell lines
                                |   *CELLCOUNT* is the number of cells tested that had nonzero signals.
                                | Provenance: wgEncodeRegTfbsClusteredV2 UCSC table
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
aa_length         STRING        The length of CDS in terms of number of amino acids
biotype           STRING        The type of transcript (e.g.protein-coding, pseudogene, rRNA etc.)
impact            STRING        Impacts due to variation (ref.impact category)
impact_severity   STRING        Severity of the impact based on the impact column value (ref.impact category)
polyphen_pred     STRING        | Impact of the SNP as given by PolyPhen (VEP only) 
                                | benign, possibly_damaging, probably_damaging, unknown
polyphen_scores   FLOAT         Polyphen score reflecting severity (higher the impact, *higher* the score)
sift_pred         STRING        | Impact of the SNP as given by SIFT (VEP only)
                                | neutral, deleterious
sift_scores       FLOAT         SIFT prob. scores reflecting severity (Higher the impact, *lower* the score)
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
HIGH              - exon_deleted
                  - frame_shift
                  - splice_acceptor
                  - splice_donor
                  - start_loss
                  - stop_gain
                  - stop_loss
                  - non_synonymous_start
MED               - non_syn_coding
                  - inframe_codon_gain
                  - inframe_codon_loss
                  - inframe_codon_change
                  - codon_change_del
                  - codon_change_ins
                  - UTR_5_del
                  - UTR_3_del
                  - other_splice_variant
                  - mature_miRNA
                  - regulatory_region
                  - TF_binding_site
                  - regulatory_region_ablation
                  - regulatory_region_amplification
                  - TFBS_ablation
                  - TFBS_amplification 
LOW               - synonymous_stop
                  - synonymous_coding
                  - UTR_5_prime
                  - UTR_3_prime
                  - intron
                  - CDS
                  - upstream
                  - downstream
                  - intergenic
                  - intragenic
                  - gene
                  - transcript
                  - exon
                  - start_gain
                  - synonymous_start
                  - intron_conserved
                  - nc_transcript
                  - NMD_transcript
                  - transcript_codon_change
                  - incomplete_terminal_codon
                  - nc_exon
                  - transcript_ablation
                  - transcript_amplification
                  - feature elongation
                  - feature truncation   
================  =======================================



