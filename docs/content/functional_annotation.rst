#############################
Annotation with snpEff or VEP
#############################

`GEMINI` depends upon external tools to predict the functional consequence of variants in a VCF file.
We currently support annotations produced by either `SnpEff <http://snpeff.sourceforge.net/>`_
or `VEP <http://www.ensembl.org/info/docs/variation/vep/index.html>`_.

.. note::
	
	Versions tested: VEP versions 73 through 75 and core SnpEff versions 3.0 through 3.6.
	GEMINI supports ENSEMBL annotations hence users are expected to download genome databases
	for these tools as represented in the examples below.
	
.. note::
 	
	Version support would be subsequently updated here, as we test along 
	and add or edit changes available with the latest version of these tools.
	

Recommended instructions for annotating existing VCF files with these tools are summarized here.

Stepwise installation and usage of VEP
=======================================

Download the Variant Effect Predictor "standalone perl script" from 
`Ensembl <http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html>`_. 
You can choose a specific version of VEP to download here

Example:

.. code-block:: bash
	
	Download version 74

Untar the tarball into the current directory

.. code-block:: bash
    
    $ tar -zxvf variant_effect_predictor.tar.gz

This will create the variant_effect_predictor directory. Now do the following for install:

.. code-block:: bash

    $ cd variant_effect_predictor
    $ perl INSTALL.pl [options]


By default this would install the API's, bioperl-1.2.3 and the cache files (in 
the $HOME/.vep directory).

Homebrew VEP installation
-------------------------

If you are a `Homebrew <http://brew.sh/>`_ or `Linuxbrew <http://brew.sh/linuxbrew/>`_
user, there is an automated recipe to install the main VEP script and plugins in
the `CloudBioLinux homebrew repository <https://github.com/chapmanb/homebrew-cbl>`_:

.. code-block:: bash

    $ brew tap chapmanb/cbl
    $ brew update
    $ brew install vep

Manual installation of VEP
--------------------------

For those (e.g mac users) who have a problem installing through this install script, try a 
manual installation of the API's, BioPerl-1.2.3 and set all pre-requisites for running VEP 
(``DBI`` and ``DBD::mysql`` modules required). The appropriate pre-build caches should be 
downloaded for Human to the ``$HOME/.vep`` directory and then untar.

You may follow instructions at http://www.ensembl.org/info/docs/api/api_installation.html
which provides alternate options for the API installation and additional tips for windows/mac 
users. It also has information for setting up your environment to run VEP.

Example download of the cache files

.. code-block:: bash

	
	$ wget ftp://ftp.ensembl.org/pub/release-73/variation/VEP/homo_sapiens_vep_73.tar.gz
	
You may change the release date in this example to get the appropriate cache files for your version
of VEP that you have installed.

Example

.. code-block:: bash
	
	
	$ wget ftp://ftp.ensembl.org/pub/release-74/variation/VEP/homo_sapiens_vep_74.tar.gz

Cache requires the ``gzip`` and ``zcat`` utilities. VEP uses ``zcat`` 
to decompress cached files. For systems where zcat may not be installed or may 
not work, the following option needs to be added along with the ``--cache`` option:

.. code-block:: bash

    --compress "gunzip -c"

Running VEP
-----------
You may now run VEP as:

.. code-block:: bash

 
    $ perl variant_effect_predictor.pl [OPTIONS]
    
.. note::

We recommend running VEP with the following options as currently we support 
VEP fields specified as below:

.. code-block:: bash

    $ perl variant_effect_predictor.pl -i example.vcf \
	--cache \
	--sift b \
	--polyphen b \
	--symbol \
	--numbers \
	--biotype \
	--total_length \
	-o output \
	--vcf \
	--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
		
    
A documentation for the above specified options may be found at 
http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html

As of GEMINI version 0.8.0, you can also run VEP with additional fields, which
will be automatically added to the variants table as columns. As an example,
run VEP on your VCF with the `dbNSFP <https://github.com/ensembl-variation/VEP_plugins/blob/master/dbNSFP.pm>`_
and `LOFTEE <https://github.com/konradjk/loftee>`_ plugins to annotate potential high
impact variations:

.. code-block:: bash

    $ variant_effect_predictor.pl --sift b --polyphen b --symbol --numbers --biotype \
    --total_length --canonical --ccds \
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,LoF,LoF_filter,LoF_flags \
    --plugin dbNSFP,/path/to/dbNSFP_v2.5.gz,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index \
    --plugin LoF,human_ancestor_fa:/path/to/human_ancestor.fa

Feeding this into GEMINI produces a variants table with columns for each of the
additional VEP metrics. The annotation loader names each column by prefixing
``vep_`` to the origin VEP name, so select on ``vep_radialsvm_score`` or
``vep_lof_filter`` in the final database.

Stepwise installation and usage of SnpEff
=======================================


.. note::

    Basic Requirements: Java v1.7 or later; at least 4GB of memory


Download the supported versions of SnpEff from http://snpeff.sourceforge.net/download.html 

Example:

.. code-block:: bash

    $ wget http://sourceforge.net/projects/snpeff/files/snpEff_v3_6_core.zip

.. note::

    SnpEff should be installed preferably in ``snpEff`` directory in your 
    home directory. Else, you must update the ``data_dir`` parameter in 
    your snpEff.config file. For e.g. if the installation of snpEff has been done 
    in ``~/src`` instead of ``~/`` then change the data_dir parameter in 
    snpEff.config to ``data_dir = ~/src/snpEff/data/``


Unzip the downloaded package.

.. code-block:: bash

    $ unzip snpEff_v3_6_core.zip


Change to the ``snpEff`` directory and download the genome database.

.. code-block:: bash

    $ cd snpEff_v3_6_core
    $ java -jar snpEff.jar download GRCh37.69


Unzip the downloaded genome database. This will create and place the genome 
in the 'data' directory                                                                                  

.. code-block:: bash

    $ unzip snpEff_v3_6_GRCh37.69.zip

                                                                                            
To annotate a vcf using snpEff, use the ``default`` options as below:


.. note::

    Memory options for the run may be specified as ``-Xmx4G`` (4GB)

.. code-block:: bash


    $ java -Xmx4G -jar snpEff.jar -i vcf -o vcf GRCh37.69 example.vcf > example_snpeff.vcf


If running from a directory different from the installation directory, the 
complete path needs to be specified as,  e.g.:

.. code-block:: bash

    $ java -Xmx4G -jar path/to/snpEff/snpEff.jar -c path/to/snpEff/snpEff.config GRCh37.69 path/to/example.vcf > example_snpeff.vcf

.. note::

	When using the latest versions of snpEff (e.g. 4.1) annotate your VCF with the additional parameters ``-classic`` and ``-formatEff``.
	This would ensure proper loading of the gene info columns in the variants table.

Columns populated by snpEff/VEP tools
=====================================

The following variant consequence columns in the variant/variant_impacts table,
are populated with these annotations, which are otherwise set to null.

* anno_id
* gene
* transcript
* exon
* is_exonic
* is_lof
* is_coding
* codon_change
* aa_change
* aa_length
* biotype
* impact
* impact_so
* impact_severity
* polyphen_pred
* polyphen_score
* sift_pred
* sift_score


Standardizing ``impact`` definitions for GEMINI
===============================================

GEMINI uses slightly modified impact terms (for ease) to describe the functional consequence of a given variant as provided by
snpEff/VEP.

The table below shows the alternate `GEMINI` terms used for `snpEff/VEP`.

=============================================       ===================================    =================================================
GEMINI terms                                        snpEff terms                           VEP terms (uses SO by default)
=============================================       ===================================    =================================================
splice_acceptor                                     SPLICE_SITE_ACCEPTOR                   splice_acceptor_variant_
splice_donor                                        SPLICE_SITE_DONOR                      splice_donor_variant_
stop_gain                                           STOP_GAINED                            stop_gained_
stop_loss                                           STOP_LOST                              stop_lost_
frame_shift                                         FRAME_SHIFT                            frameshift_variant_
start_loss                                          START_LOST                             null
exon_deleted                                        EXON_DELETED                           null
non_synonymous_start                                NON_SYNONYMOUS_START                   null
transcript_codon_change                             null                                   initiator_codon_variant_
chrom_large_del                                     CHROMOSOME_LARGE_DELETION              null
rare_amino_acid                                     RARE_AMINO_ACID                        null
non_syn_coding                                      NON_SYNONYMOUS_CODING                  missense_variant_
inframe_codon_gain                                  CODON_INSERTION                        inframe_insertion_
inframe_codon_loss                                  CODON_DELETION                         inframe_deletion_
inframe_codon_change                                CODON_CHANGE                           null
codon_change_del                                    CODON_CHANGE_PLUS_CODON_DELETION       null
codon_change_ins                                    CODON_CHANGE_PLUS_CODON_INSERTION      null
UTR_5_del                                           UTR_5_DELETED                          null
UTR_3_del                                           UTR_3_DELETED                          null
splice_region                                       SPLICE_SITE_REGION                     splice_region_variant_
mature_miRNA                                        null                                   mature_miRNA_variant_
regulatory_region                                   null                                   regulatory_region_variant_
TF_binding_site                                     null                                   TF_binding_site_variant_
regulatory_region_ablation                          null                                   regulatory_region_ablation_
regulatory_region_amplification                     null                                   regulatory_region_amplification_
TFBS_ablation                                       null                                   TFBS_ablation_
TFBS_amplification                                  null                                   TFBS_amplification_
synonymous_stop                                     SYNONYMOUS_STOP                        stop_retained_variant_
synonymous_coding                                   SYNONYMOUS_CODING                      synonymous_variant_
UTR_5_prime                                         UTR_5_PRIME                            5_prime_UTR_variant_
UTR_3_prime                                         UTR_3_PRIME                            3_prime_UTR_variant_
intron                                              INTRON                                 intron_variant_
CDS                                                 CDS                                    coding_sequence_variant_
upstream                                            UPSTREAM                               upstream_gene_variant_
downstream                                          DOWNSTREAM                             downstream_gene_variant_
intergenic                                          INTERGENIC                             intergenic_variant_
intergenic_conserved                                INTERGENIC_CONSERVED                   null
intragenic                                          INTRAGENIC                             null
gene                                                GENE                                   null
transcript                                          TRANSCRIPT                             null
exon                                                EXON                                   null
start_gain                                          START_GAINED                           null
synonymous_start                                    SYNONYMOUS_START                       null
intron_conserved                                    INTRON_CONSERVED                       null
nc_transcript                                       null                                   nc_transcript_variant_ (should have been returned by VEP as: non_coding_transcript_variant)
NMD_transcript                                      null                                   NMD_transcript_variant_
incomplete_terminal_codon                           null                                   incomplete_terminal_codon_variant_
nc_exon                                             null                                   non_coding_exon_variant_ (should have been returned by VEP as: non_coding_transcript_exon_variant)
transcript_ablation                                 null                                   transcript_ablation_
transcript_amplification                            null                                   transcript_amplification_
feature elongation                                  null                                   feature_elongation_
feature truncation                                  null                                   feature_truncation_
=============================================       ===================================    =================================================
*Note: "null" refers to the absence of the corresponding term in the alternate database*



SO ``impact`` definitions in GEMINI
===============================================

The below table shows the `Sequence Ontology (SO)` term mappings for the GEMINI impacts, which is otherwise contained in the
``impact_so`` column of the variants/variant_impacts table of the GEMINI database. The last column shows the severity terms
defined in GEMINI for these impacts.

=============================================     =======================================================     =================
GEMINI terms (column: impact)                     Sequence Ontology terms (column: impact_so)                 Impact severity
=============================================     =======================================================     =================
splice_acceptor                                   splice_acceptor_variant_                                    HIGH
splice_donor                                      splice_donor_variant_                                       HIGH
stop_gain                                         stop_gained_                                                HIGH
stop_loss                                         stop_lost_                                                  HIGH
frame_shift                                       frameshift_variant_                                         HIGH
start_loss                                        start_lost_                                                 HIGH
exon_deleted                                      exon_loss_variant_                                          HIGH
non_synonymous_start                              initiator_codon_variant_                                    HIGH
transcript_codon_change                           initiator_codon_variant_                                    HIGH
chrom_large_del                                   chromosomal_deletion_                                       HIGH
rare_amino_acid                                   rare_amino_acid_variant_                                    HIGH
non_syn_coding                                    missense_variant_                                           MED
inframe_codon_gain                                inframe_insertion_                                          MED
inframe_codon_loss                                inframe_deletion_                                           MED
inframe_codon_change                              coding_sequence_variant_                                    MED
codon_change_del                                  disruptive_inframe_deletion_                                MED
codon_change_ins                                  disruptive_inframe_insertion_                               MED
UTR_5_del                                         5_prime_UTR_truncation_ + exon_loss_variant_                MED
UTR_3_del                                         3_prime_UTR_truncation_ + exon_loss_variant_                MED
splice_region                                     splice_region_variant_                                      MED
mature_miRNA                                      mature_miRNA_variant_                                       MED
regulatory_region                                 regulatory_region_variant_                                  MED
TF_binding_site                                   TF_binding_site_variant_                                    MED
regulatory_region_ablation                        regulatory_region_ablation_                                 MED
regulatory_region_amplification                   regulatory_region_amplification_                            MED
TFBS_ablation                                     TFBS_ablation_                                              MED
TFBS_amplification                                TFBS_amplification_                                         MED
synonymous_stop                                   stop_retained_variant_                                      LOW
synonymous_coding                                 synonymous_variant_                                         LOW
UTR_5_prime                                       5_prime_UTR_variant_                                        LOW
UTR_3_prime                                       3_prime_UTR_variant_                                        LOW
intron                                            intron_variant_                                             LOW
CDS                                               coding_sequence_variant_                                    LOW
upstream                                          upstream_gene_variant_                                      LOW
downstream                                        downstream_gene_variant_                                    LOW
intergenic                                        intergenic_variant_                                         LOW
intergenic_conserved                              conserved_intergenic_variant_                               LOW
intragenic                                        intragenic_variant_                                         LOW
gene                                              gene_variant_                                               LOW
transcript                                        transcript_variant_                                         LOW
exon                                              exon_variant_                                               LOW
start_gain                                        5_prime_UTR_premature_start_codon_gain_variant_             LOW
synonymous_start                                  start_retained_variant_                                     LOW
intron_conserved                                  conserved_intron_variant_                                   LOW
nc_transcript                                     nc_transcript_variant_                                      LOW
NMD_transcript                                    NMD_transcript_variant_                                     LOW
incomplete_terminal_codon                         incomplete_terminal_codon_variant_                          LOW
nc_exon                                           non_coding_exon_variant_                                    LOW
transcript_ablation                               transcript_ablation_                                        LOW
transcript_amplification                          transcript_amplification_                                   LOW
feature elongation                                feature_elongation_                                         LOW
feature truncation                                feature_truncation_                                         LOW
=============================================     =======================================================     =================

.. _3_prime_UTR_truncation: http://sequenceontology.org/browser/current_svn/term/SO:0002015
.. _3_prime_UTR_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001624
.. _5_prime_UTR_premature_start_codon_gain_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001988
.. _5_prime_UTR_truncation: http://sequenceontology.org/browser/current_svn/term/SO:0002013
.. _5_prime_UTR_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001623
.. _chromosomal_deletion: http://sequenceontology.org/browser/current_svn/term/SO:1000029
.. _coding_sequence_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001580
.. _conserved_intergenic_variant: http://sequenceontology.org/browser/current_svn/term/SO:0002017
.. _conserved_intron_variant: http://sequenceontology.org/browser/current_svn/term/SO:0002018
.. _disruptive_inframe_deletion: http://sequenceontology.org/browser/current_svn/term/SO:0001826
.. _disruptive_inframe_insertion: http://sequenceontology.org/browser/current_svn/term/SO:0001824
.. _downstream_gene_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001632
.. _exon_loss_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001572
.. _exon_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001791
.. _feature_elongation: http://sequenceontology.org/browser/current_svn/term/SO:0001907
.. _feature_truncation: http://sequenceontology.org/browser/current_svn/term/SO:0001906
.. _frameshift_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001589
.. _gene_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001564
.. _incomplete_terminal_codon_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001626
.. _inframe_deletion: http://sequenceontology.org/browser/current_svn/term/SO:0001822
.. _inframe_insertion: http://sequenceontology.org/browser/current_svn/term/SO:0001821
.. _initiator_codon_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001582
.. _intergenic_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001628
.. _intragenic_variant: http://sequenceontology.org/browser/current_svn/term/SO:0002011
.. _intron_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001627
.. _mature_miRNA_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001620
.. _missense_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001583
.. _NMD_transcript_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001621
.. _non_coding_exon_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001792
.. _rare_amino_acid_variant: http://sequenceontology.org/browser/current_svn/term/SO:0002008
.. _regulatory_region_ablation: http://sequenceontology.org/browser/current_svn/term/SO:0001894
.. _regulatory_region_amplification: http://sequenceontology.org/browser/current_svn/term/SO:0001891
.. _regulatory_region_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001566
.. _splice_acceptor_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001574
.. _splice_donor_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001575
.. _splice_region_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001630
.. _start_lost: http://sequenceontology.org/browser/current_svn/term/SO:0002012
.. _start_retained_variant: http://sequenceontology.org/browser/current_svn/term/SO:0002019
.. _stop_gained: http://sequenceontology.org/browser/current_svn/term/SO:0001587
.. _stop_lost: http://sequenceontology.org/browser/current_svn/term/SO:0001578
.. _stop_retained_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001567
.. _synonymous_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001819
.. _TF_binding_site_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001782
.. _TFBS_ablation: http://sequenceontology.org/browser/current_svn/term/SO:0001895
.. _TFBS_amplification: http://sequenceontology.org/browser/current_svn/term/SO:0001892
.. _transcript_ablation: http://sequenceontology.org/browser/current_svn/term/SO:0001893
.. _transcript_amplification: http://sequenceontology.org/browser/current_svn/term/SO:0001889
.. _transcript_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001576
.. _upstream_gene_variant: http://sequenceontology.org/browser/current_svn/term/SO:0001631
.. _nc_transcript_variant: http://www.sequenceontology.org/browser/current_svn/term/SO:0001619
