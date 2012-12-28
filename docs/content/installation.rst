############
Installation
############


Software installation
=====================
To install you should download the latest source code from GitHub, either by going to::

    http://github.com/arq5x/gemini

and clicking on "Downloads", or by cloning the git repository with::

    $ git clone https://github.com/arq5x/gemini.git

Once you have the source code, run::

    $ cd gemini
    $ sudo python setup.py install

to install it. If you don't have permission to install it in the default 
directory, you can simply build the source in-place and use the package 
from the git repository::

    $  python setup.py build_ext --inplace


Installing annotation files
===========================
One of the more appealing features in ``gemini`` is that it automatically 
annotates variants in a VCF file with several genome annotations.  However, 
you must first install these data files on your system. It's easy enough --- 
you just need to run the following script and tell it in which what full path 
you'd like to install the necessary data files. The recommended path is 
``/usr/local/share``, but you can install the data files wherever you want.

    gemini/install-data.py /usr/local/share/



---------------------------
Functional annotation tools
---------------------------
`gemini` depends upon external tools to predict the functional consequence of variants in a VCF file.
We currently support annotations produced by both `SnpEff <http://snpeff.sourceforge.net/>`_ 
and `VEP <http://useast.ensembl.org/info/docs/variation/vep/index.html>`_.  
Recommended instructions for annotating existing VCF files with these tools are available here.  
In addition, we have attempted to standardize the terms used to describe the functional consequence of a given variant, 
as each annotation tool uses different vocabulary.


The variant consequence columns in the variant table are populated either by `snpEff` or `VEP` as defined by the user using the `-t` option while running pop load 
(To populate these columns the input VCF file should have been annotated either by `snpEff` or `VEP`)::

	$ pop load -v my.vcf -t VEP -d my.db
	$ pop load -v my.vcf -t snpEFF -d my.db

By default the following columns in the variant table would be set to null:

* anno_id
* gene
* affected_gene
* affected_transcript
* affected_exon
* is_exonic
* is_lof
* is_coding
* codon_change
* aa_change
* aa_length
* biotype
* most_severe_impact
* impact_severity
* polyphen_pred
* polyphen_score
* sift_pred
* sift_score

Impacts
.......
The table below shows the alternate `gemini` terms for the consequences from `snpEff` and `VEP`, for SQL queries. 
The last column represents the severity terms associated with the impacts:

=============================================       ===================================    =====================================================     ================
Gemini terms                                        snpEff terms                           VEP terms                                                 Impact severity
=============================================       ===================================    =====================================================     ================
splice_acceptor                                     SPLICE_SITE_ACCEPTOR                   splice_acceptor_variant                                   HIGH
splice_donor                                        SPLICE_SITE_DONOR                      splice_donor_variant                                      HIGH
stop_gain                                           STOP_GAINED                            stop_gained                                               HIGH
stop_loss                                           STOP_LOST                              stop_lost                                                 HIGH
frame_shift                                         FRAME_SHIFT                            frameshift_variant                                        HIGH
start_loss                                          START_LOST                             null                                                      HIGH
exon_deleted                                        EXON_DELETED                           null                                                      HIGH
non_synonymous_start                                NON_SYNONYMOUS_START                   null                                                      HIGH
non_syn_coding                                      NON_SYNONYMOUS_CODING                  missense_variant                                          MED
inframe_codon_gain                                  CODON_INSERTION                        inframe_insertion                                         MED
inframe_codon_loss                                  CODON_DELETION                         inframe_deletion                                          MED
inframe_codon_change                                CODON_CHANGE                           null                                                      MED
codon_change_del                                    CODON_CHANGE_PLUS_CODON_DELETION       null                                                      MED
codon_change_ins                                    CODON_CHANGE_PLUS_CODON_INSERTION      null                                                      MED
UTR_5_del                                           UTR_5_DELETED                          null                                                      MED
UTR_3_del                                           UTR_3_DELETED                          null                                                      MED
other_splice_variant                                null                                   splice_region_variant                                     MED
mature_miRNA                                        null                                   mature_miRNA_variant                                      MED
regulatory_region                                   null                                   regulatory_region_variant                                 MED
TF_binding_site                                     null                                   TF_binding_site_variant                                   MED
regulatory_region_ablation                          null                                   regulatory_region_ablation                                MED
regulatory_region_amplification                     null                                   regulatory_region_amplification                           MED
TFBS_ablation                                       null                                   TFBS_ablation                                             MED
TFBS_amplification                                  null                                   TFBS_amplification                                        MED                                         
synonymous_stop                                     SYNONYMOUS_STOP                        stop_retained_variant                                     LOW
synonymous_coding                                   SYNONYMOUS_CODING                      synonymous_variant                                        LOW
UTR_5_prime                                         UTR_5_PRIME                            5_prime_UTR_variant                                       LOW
UTR_3_prime                                         UTR_3_PRIME                            3_prime_UTR_variant                                       LOW
intron                                              INTRON                                 intron_variant                                            LOW
CDS                                                 CDS                                    coding_sequence_variant                                   LOW
upstream                                            UPSTREAM                               upstream_gene_variant					                 LOW              
downstream                                          DOWNSTREAM                             downstream_gene_variant                                   LOW
intergenic                                          INTERGENIC, INTERGENIC_CONSERVED       intergenic_variant                                        LOW
intragenic                                          INTRAGENIC                             null                                                      LOW
gene                                                GENE                                   null                                                      LOW
transcript                                          TRANSCRIPT                             null                                                      LOW   
exon                                                EXON                                   null                                                      LOW
start_gain                                          START_GAINED                           null                                                      LOW
synonymous_start                                    SYNONYMOUS_START                       null                                                      LOW
intron_conserved                                    INTRON_CONSERVED                       null                                                      LOW
nc_transcript                                       null                                   nc_transcript_variant                                     LOW
NMD_transcript                                      null                                   NMD_transcript_variant                                    LOW
transcript_codon_change                             null                                   initiator_codon_variant                                   LOW
incomplete_terminal_codon                           null                                   incomplete_terminal_codon_variant                         LOW
nc_exon                                             null                                   non_coding_exon_variant                                   LOW
transcript_ablation                                 null                                   transcript_ablation                                       LOW
transcript_amplification                            null                                   transcript_amplification                                  LOW
feature elongation                                  null                                   feature elongation                                        LOW
feature truncation                                  null                                   feature truncation                                        LOW
=============================================       ===================================    =====================================================     ================
*Note: "null" refers to the absence of the corresponding term in the alternate database* 