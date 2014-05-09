##########################
The Gemini database schema
##########################


The ``variants`` table
----------------------


Core VCF fields
...............
   
========================  ========      ==============================================================================================
column_name               type          notes
========================  ========      ==============================================================================================
chrom                     STRING        The chromosome on which the variant resides (from VCF ``CHROM`` field).
start                     INTEGER       The 0-based start position. (from VCF ``POS`` field, but converted to 0-based coordinates)
end                       INTEGER       The 1-based end position. (from VCF ``POS`` field, yet inferred based on the size of the variant)
vcf_id                    STRING        The VCF ``ID`` field.
variant_id                INTEGER       PRIMARY_KEY
anno_id                   INTEGER       Variant transcript number for the most severely affected transcript
ref                       STRING        Reference allele (from VCF ``REF`` field)
alt                       STRING        Alternate alele for the variant (from VCF ``ALT`` field)
qual                      INTEGER       Quality score for the assertion made in ALT (from VCF ``QUAL`` field)
filter                    STRING        A string of filters passed/failed in variant calling (from VCF ``FILTER`` field)
========================  ========      ==============================================================================================



Variant and PopGen info
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
type                      STRING        | The type of variant.
                                        | Any of: [*snp*, *indel*]
sub_type                  STRING        | The variant sub-type.
                                        | If ``type`` is *snp*:   [*ts*, (transition), *tv* (transversion)]
                                        | If ``type`` is *indel*: [*ins*, (insertion), *del* (deletion)]
call_rate                 FLOAT         The fraction of samples with a valid genotype
num_hom_ref               INTEGER       The total number of of homozygotes for the reference (``ref``) allele
num_het                   INTEGER       The total number of heterozygotes observed.
num_hom_alt               INTEGER       The total number of homozygotes for the reference (``alt``) allele
num_unknown               INTEGER       The total number of of unknown genotypes
aaf                       FLOAT         The observed allele frequency for the alternate allele
hwe                       FLOAT         The Chi-square probability of deviation from HWE (assumes random mating)
inbreeding_coeff          FLOAT         The inbreeding co-efficient that expresses the likelihood of effects due to inbreeding
pi                        FLOAT         The computed nucleotide diversity (pi) for the site
========================  ========      ==============================================================================================



Genotype information
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
gts                       BLOB          | A compressed binary vector of sample genotypes (e.g., "A/A", "A|G", "G/G")
                                        | - Extracted from the VCF ``GT`` genotype tag.
gt_types                  BLOB          | A compressed binary vector of numeric genotype "types" (e.g., 0, 1, 2)
                                        | - Inferred from the VCF ``GT`` genotype tag.
gt_phases                 BLOB          | A compressed binary vector of sample genotype phases (e.g., False, True, False)
                                        | - Extracted from the VCF ``GT`` genotype tag's allele delimiter
                                        |   e.g., ``A/G`` means an unphased genotype. Value is **FALSE**.
                                        |   e.g., ``A|G`` means a phased genotype. Value is **TRUE**.
gt_depths                 BLOB          | A compressed binary vector of the depth of aligned sequence observed for each sample
                                        | - Extracted from the VCF ``DP`` genotype tag.
gt_ref_depths             BLOB          | A compressed binary vector of the depth of reference alleles observed for each sample
                                        | - Extracted from the VCF ``AD`` genotype tag.
gt_alt_depths             BLOB          | A compressed binary vector of the depth of alternate alleles observed for each sample
                                        | - Extracted from the VCF ``AD`` genotype tag.
gt_quals                  BLOB          | A compressed binary vector of the genotype quality (PHRED scale) estimates for each sample
                                        | - Extracted from the VCF ``GQ`` genotype tag.
========================  ========      ==============================================================================================



Gene information
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
gene                      STRING        Corresponding gene name of the highly affected transcript
transcript                STRING        | The variant transcript that was most severely affected
                                        | (for two equally affected transcripts,the protein_coding
										biotype is prioritized (snpEff/VEP)
is_exonic                 BOOL          Does the variant affect an exon for >= 1transcript?
is_coding                 BOOL          Does the variant fall in a coding region (excl. 3' & 5' UTRs) for >= 1 transcript?
is_lof                    BOOL          Based on the value of the impact col, is the variant LOF for >= transcript?
exon                      STRING        Exon information for the severely affected transcript
codon_change              STRING        What is the codon change?
aa_change                 STRING        What is the amino acid change (for an snp)?
aa_length                 STRING        The length of CDS in terms of number of amino acids (``only snpEff``)
biotype                   STRING        The 'type' of the severely affected transcript (e.g.protein-coding, pseudogene, rRNA etc.) (``only snpEff``)
impact                    STRING        The consequence of the most severely affected transcript
impact_so                 STRING        The Sequence ontology term for the most severe consequence
impact_severity           STRING        Severity of the highest order observed for the variant
polyphen_pred             STRING        Polyphen predictions for the snps for the severely affected transcript (``only VEP``) 
polyphen_score            FLOAT         Polyphen scores for the severely affected transcript (``only VEP``)
sift_pred                 STRING        SIFT predictions for the snp's for the most severely affected transcript (``only VEP``)
sift_score                FLOAT         SIFT scores for the predictions (``only VEP``)
pfam_domain               STRING        Pfam protein domain that the variant affects
========================  ========      ==============================================================================================



Optional VCF INFO fields
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
anc_allele                STRING        The reported ancestral allele if there is one.
rms_bq                    FLOAT         The RMS base quality at this position.
cigar                     STRING        CIGAR string describing how to align an alternate allele to the reference allele.
depth                     INTEGER       The number of aligned sequence reads that led to this variant call
strand_bias               FLOAT         Strand bias at the variant position
rms_map_qual              FLOAT         RMS mapping quality, a measure of variance of quality scores
in_hom_run                INTEGER       Homopolymer runs for the variant allele
num_mapq_zero             INTEGER       Total counts of reads with mapping quality equal to zero
num_alleles               INTEGER       Total number of alleles in called genotypes
num_reads_w_dels          FLOAT         Fraction of reads with spanning deletions
haplotype_score           FLOAT         Consistency of the site with two segregating haplotypes
qual_depth                FLOAT         Variant confidence or quality by depth
allele_count              INTEGER       Allele counts in genotypes
allele_bal                FLOAT         Allele balance for hets
info                      BLOB          Stores the ``INFO`` field of the VCF
========================  ========      ==============================================================================================



Population information
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
in_dbsnp                  BOOL          | Is this variant found in dbSnp (build 137)?
                                        | 0 : Absence of the variant in dbsnp
                                        | 1 : Presence of the variant in dbsnp
rs_ids                    STRING        | A comma-separated list of rs ids for variants present in dbsnp
in_hm2                    BOOL          Whether the variant was part of HapMap2.
in_hm3                    BOOL          Whether the variant was part of HapMap3.
in_esp                    BOOL          Presence/absence of the variant in the ESP project data
in_1kg                    BOOL          Presence/absence of the variant in the 1000 genome project data
aaf_esp_ea                FLOAT         Minor Allele Frequency of the variant for European Americans in the ESP project
aaf_esp_aa                FLOAT         Minor Allele Frequency of the variant for African Americans in the ESP project
aaf_esp_all               FLOAT         Minor Allele Frequency of the variant w.r.t both groups in the ESP project
aaf_1kg_amr               FLOAT         Allele Frequency of the variant for samples in AMR based on AC/AN (1000g project)
aaf_1kg_asn               FLOAT         Allele frequency of the variant for samples in ASN based on AC/AN (1000g project)
aaf_1kg_afr               FLOAT         Allele frequency of the variant for samples in AFR based on AC/AN (1000g project)
aaf_1kg_eur               FLOAT         Allele Frequency of the variant for samples in EUR based on AC/AN (1000g project)
aaf_1kg_all               FLOAT         Global allele frequency (based on AC/AN) (1000g project) 
========================  ========      ==============================================================================================



Disease phenotype info (from ClinVar).
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
in_omim                   BOOL          | 0 : Absence of the variant in OMIM database
                                        | 1 : Presence of the variant in OMIM database
clinvar_sig               STRING        | The clinical significance scores for each
                                        | of the variant according to ClinVar:
                                        | *unknown*, *untested*, *non-pathogenic*
                                        | *probable-non-pathogenic*, *probable-pathogenic*
                                        | *pathogenic*, *drug-response*, *histocompatibility*
                                        | *other*
clinvar_disease_name      STRING        The name of the disease to which the variant is relevant
clinvar_dbsource          STRING        Variant Clinical Channel IDs
clinvar_dbsource_id       STRING        The record id in the above database
clinvar_origin            STRING        | The type of variant.
                                        | Any of:
                                        | *unknown*, *germline*, *somatic*,
                                        | *inherited*, *paternal*, *maternal*,
                                        | *de-novo*, *biparental*, *uniparental*, 
                                        | *not-tested*, *tested-inconclusive*, 
                                        | *other*
clinvar_dsdb              STRING        Variant disease database name
clinvar_dsdbid            STRING        Variant disease database ID
clinvar_disease_acc       STRING        Variant Accession and Versions
clinvar_in_locus_spec_db  BOOL          Submitted from a locus-specific database?
clinvar_on_diag_assay     BOOL          Variation is interrogated in a clinical diagnostic assay?
========================  ========      ==============================================================================================



Genome annotations
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
exome_chip                BOOL          Whether an SNP is on the Illumina HumanExome Chip
cyto_band                 STRING        Chromosomal cytobands that a variant overlaps
rmsk                      STRING        | A comma-separated list of RepeatMasker annotations that the variant overlaps.
                                        | Each hit is of the form: ``name_class_family``
in_cpg_island             BOOL          | Does the variant overlap a CpG island?.
                                        | Based on UCSC: Regulation > CpG Islands > cpgIslandExt 
in_segdup                 BOOL          | Does the variant overlap a segmental duplication?.
                                        | Based on UCSC: Variation&Repeats > Segmental Dups > genomicSuperDups track
is_conserved              BOOL          | Does the variant overlap a conserved region?
                                        | Based on the 29-way mammalian conservation study
gerp_bp_score             FLOAT         | GERP conservation score. 
                                        | Only populated if the ``--load-gerp-bp`` option is used when loading. 
                                        | Higher scores reflect greater conservation. **At base-pair resolution**.
                                        | Details: http://mendel.stanford.edu/SidowLab/downloads/gerp/
gerp_element_pval         FLOAT         | GERP elements P-val 
                                        | Lower P-values scores reflect greater conservation. **Not at base-pair resolution**.  
                                        | Details: http://mendel.stanford.edu/SidowLab/downloads/gerp/
recomb_rate               FLOAT         | Returns the mean recombination rate at the variant site
                                        | Based on HapMapII_GRCh37 genetic map
cadd_raw                  FLOAT         | Raw ``CADD`` scores for scoring deleteriousness of SNV's in the human genome
                                        | Details: http://www.ncbi.nlm.nih.gov/pubmed/24487276
cadd_scaled               FLOAT         | Scaled ``CADD`` scores (Phred like) for scoring deleteriousness of SNV's
                                        | Details: http://www.ncbi.nlm.nih.gov/pubmed/24487276
========================  ========      ==============================================================================================

**Note:**
``CADD`` scores (http://cadd.gs.washington.edu/) are Copyright 2013 University of Washington and Hudson-Alpha Institute for Biotechnology 
(all rights reserved) but are freely available for all academic, non-commercial applications. For commercial licensing information contact 
Jennifer McCullar (mccullaj@uw.edu).



Variant error assessment
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
grc                       STRING        | Association with patch and fix regions from the Genome Reference Consortium:
                                        | http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
                                        | Identifies potential problem regions associated with variant calls.
                                        | Built with `annotation_provenance/make-ncbi-grc-patches.py`
gms_illumina              FLOAT         | Genome Mappability Scores (GMS) for Illumina error models
                                        | Provides low GMS scores (< 25.0 in any technology) from:
                                        | http://sourceforge.net/apps/mediawiki/gma-bio/index.php?title=Download_GMS
                                        | #Download_GMS_by_Chromosome_and_Sequencing_Technology
                                        | Input VCF for annotations prepared with:
                                        | https://github.com/chapmanb/bcbio.variation/blob/master/src/bcbio/variation/utils/gms.clj
gms_solid                 FLOAT         Genome Mappability Scores with SOLiD error models
gms_iontorrent            FLOAT         Genome Mappability Scores with IonTorrent error models
in_cse                    BOOL          | Is a variant in an error prone genomic position,
                                        | using CSE: Context-Specific Sequencing Errors 
                                        | https://code.google.com/p/discovering-cse/
                                        | http://www.biomedcentral.com/1471-2105/14/S5/S1
========================  ========      ==============================================================================================



ENCODE information
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
encode_tfbs               STRING        | Comma-separated list of transcription factors that were
                                        | observed by ENCODE to bind DNA in this region.  Each hit in the list is constructed
                                        | as TF_CELLCOUNT, where:
                                        |   *TF* is the transcription factor name
                                        |   *CELLCOUNT* is the number of cells tested that had nonzero signals.
                                        | Provenance: wgEncodeRegTfbsClusteredV2 UCSC table
encode_dnaseI_cell_count  INTEGER       | Count of cell types that were observed to have DnaseI hypersensitivity.
encode_dnaseI_cell_list   STRING        | Comma separated list of cell types that were observed to have DnaseI hypersensitivity.
                                        | Provenance: Thurman, et al, *Nature*, 489, pp. 75-82, 5 Sep. 2012
encode_consensus_gm12878  STRING        | ENCODE consensus segmentation prediction for GM12878.
                                        |
                                        | CTCF: CTCF-enriched element
                                        | E:    Predicted enhancer
                                        | PF:   Predicted promoter flanking region
                                        | R:    Predicted repressed or low-activity region
                                        | TSS:  Predicted promoter region including TSS
                                        | T:    Predicted transcribed region
                                        | WE:   Predicted weak enhancer or open chromatin cis-regulatory element
						  			    | unknown: This region of the genome had no functional prediction.
encode_consensus_h1hesc   STRING        ENCODE consensus segmentation prediction for h1HESC.  See encode_consseg_gm12878 for details.       
encode_consensus_helas3   STRING        ENCODE consensus segmentation prediction for Helas3.  See encode_consseg_gm12878 for details.         
encode_consensus_hepg2    STRING        ENCODE consensus segmentation prediction for HEPG2.   See encode_consseg_gm12878 for details.          
encode_consensus_huvec    STRING        ENCODE consensus segmentation prediction for HuVEC.   See encode_consseg_gm12878 for details.        
encode_consensus_k562     STRING        ENCODE consensus segmentation prediction for k562.    See encode_consseg_gm12878 for details.
vista_enhancers           STRING        Experimentally validated human enhancers from VISTA (http://enhancer.lbl.gov/frnt_page_n.shtml)
========================  ========      ==============================================================================================



Cancer related columns
........................
========================  ========      ==============================================================================================
========================  ========      ==============================================================================================
is_somatic                BOOL          | Whether the variant is somatically acquired.
cosmic_ids                STRING        | A list of known COSMIC ids for this variant.
                                        | See: http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/
========================  ========      ==============================================================================================


|

The ``variant_impacts`` table
-----------------------------
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
aa_length         STRING        The length of CDS in terms of number of amino acids (``snpEff only``)
biotype           STRING        The type of transcript (e.g.protein-coding, pseudogene, rRNA etc.) (``SnpEff only``)
impact            STRING        Impacts due to variation (ref.impact category)
impact_so         STRING        The sequence ontology term for the impact
impact_severity   STRING        Severity of the impact based on the impact column value (ref.impact category)
polyphen_pred     STRING        | Impact of the SNP as given by PolyPhen (``VEP only``) 
                                | benign, possibly_damaging, probably_damaging, unknown
polyphen_scores   FLOAT         Polyphen score reflecting severity (higher the impact, *higher* the score) (``VEP only``)
sift_pred         STRING        | Impact of the SNP as given by SIFT (``VEP only``)
                                | neutral, deleterious
sift_scores       FLOAT         SIFT prob. scores reflecting severity (Higher the impact, *lower* the score) (``VEP only``)
================  ========      ===============================================================================

|

Details of the ``impact`` and ``impact_severity`` columns
---------------------------------------------------------

=====================  =======================================
impact severity        impacts
=====================  =======================================
HIGH                   - exon_deleted
                       - frame_shift
                       - splice_acceptor
                       - splice_donor
                       - start_loss
                       - stop_gain
                       - stop_loss
                       - non_synonymous_start
                       - transcript_codon_change
                       - rare_amino_acid
                       - chrom_large_del
MED                    - non_syn_coding
                       - inframe_codon_gain
                       - inframe_codon_loss
                       - inframe_codon_change
                       - codon_change_del
                       - codon_change_ins
                       - UTR_5_del
                       - UTR_3_del
                       - splice_region
                       - mature_miRNA
                       - regulatory_region
                       - TF_binding_site
                       - regulatory_region_ablation
                       - regulatory_region_amplification
                       - TFBS_ablation
                       - TFBS_amplification 
LOW                    - synonymous_stop
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
=====================  =======================================

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


The ``resources`` table
-----------------------

Establishes provenance of annotation resources used to create a Gemini database.

=============  ==========  ==================================================
column name    type        notes
=============  ==========  ==================================================
name           STRING      Name of the annotation type
resource       STRING      Filename of the resource, with version information
=============  ==========  ==================================================


The ``version`` table
-----------------------

Establishes which version of ``gemini`` was used to create a database.

=============  ==========  ==================================================
column name    type        notes
=============  ==========  ==================================================
version        STRING      What version of gemini was used to create the DB.
=============  ==========  ==================================================



The ``gene_detailed`` table
---------------------------

Built on version 75 of Ensembl genes

==================  ========      ===============================================================================
column_name         type          notes
==================  ========      ===============================================================================
uid                 INTEGER       PRIMARY_KEY (unique identifier for each entry in the table)
chrom               STRING        The chromosome on which the gene resides     
gene                STRING        The gene name
is_hgnc             BOOL          Flag for gene column: 0 for non HGNC symbol and 1 for HGNC symbol = TRUE
ensembl_gene_id     STRING        The ensembl gene id for the gene
transcript          STRING        The ensembl transcript id for the gene
biotype             STRING        The biotype (e.g protein coding) of the transcript
transcript_status   STRING        The status of the transcript (e.g. KNOWN, PUTATIVE etc.)
ccds_id             STRING        The consensus coding sequence transcript identifier
hgnc_id             STRING        The HGNC identifier for the gene if HGNC symbol is TRUE
entrez_id           STRING        The entrez gene identifier for the gene
cds_length          STRING        The length of CDS in bases
protein_length      STRING        The length of the transcript as the number of amino acids
transcript_start    STRING        The start position of the transcript in bases
transcript_end      STRING        The end position of the transcript in bases
strand              STRING        The strand of DNA where the gene resides
synonym             STRING        Other gene names (previous or synonyms) for the gene
rvis_pct            FLOAT         The RVIS percentile values for the gene
mam_phenotype_id    STRING        | High level mammalian phenotype ID applied to mouse phenotype descriptions
                                  | in the MGI database at http://www.informatics.jax.org/. Data taken from
								  ftp://ftp.informatics.jax.org/pub/reports/HMD_HumanPhenotype.rpt
==================  ========      ===============================================================================


The ``gene_summary`` table
---------------------------

Built on version 75 of Ensembl genes

======================  ========      ===============================================================================
column_name             type          notes
======================  ========      ===============================================================================
uid                     INTEGER       PRIMARY_KEY (unique identifier for each entry in the table)
chrom                   STRING        The chromosome on which the gene resides     
gene                    STRING        The gene name
is_hgnc                 BOOL          Flag for gene column: 0 for non HGNC symbol and 1 for HGNC symbol = TRUE
ensembl_gene_id         STRING        The ensembl gene id for the gene
hgnc_id                 STRING        The HGNC identifier for the gene if HGNC symbol is TRUE
transcript_min_start    STRING        The minimum start position of all transcripts for the gene
transcript_max_end      STRING        The maximum end position of all transcripts for the gene
strand                  STRING        The strand of DNA where the gene resides
synonym                 STRING        Other gene names (previous or synonyms) for the gene
rvis_pct                FLOAT         The RVIS percentile values for the gene
mam_phenotype_id        STRING        | High level mammalian phenotype ID applied to mouse phenotype descriptions
                                      | in the MGI database at http://www.informatics.jax.org/. Data taken from
									  ftp://ftp.informatics.jax.org/pub/reports/HMD_HumanPhenotype.rpt
in_cosmic_census        BOOL          Are mutations in the gene implicated in cancer by the cancer gene census?
======================  ========      ===============================================================================

