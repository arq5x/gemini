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

.. note::

    By default, the ``comp_hets`` tool requires phased genotypes.  If you want
    to ignore phasing in search of _putative_ compound heterozygotes, please
    see the ``--ignore-phasing`` option below.

Example usage with default parameters:


.. note::

    Each pair of consecutive lines in the output represent the two variants
    for a compound heterozygote in a give sample.  The third column,
    `comp_het_id`, tracks the distinct compound heterozygote variant pairs.

.. code-block:: bash

    $ gemini comp_hets my.db

    family  sample  comp_het_id     chrom   start   end     variant_id      anno_id ref     alt     qual    filter  type    sub_type        call_rate       in_dbsnp        rs_ids  in_omim clinvar_sig     clinvar_disease_name    clinvar_dbsource        clinvar_dbsource_id     clinvar_origin  clinvar_dsdb    clinvar_dsdbid  clinvar_disease_acc     clinvar_in_locus_spec_db        clinvar_on_diag_assay   pfam_domain     cyto_band       rmsk    in_cpg_island   in_segdup       is_conserved    gerp_bp_score   gerp_element_pval       num_hom_ref     num_het num_hom_alt     num_unknown     aaf     hwe     inbreeding_coeff        pi      recomb_rate     gene    transcript      is_exonic       is_coding       is_lof  exon    codon_change    aa_change       aa_length       biotype impact  impact_severity polyphen_pred   polyphen_score  sift_pred       sift_score      anc_allele      rms_bq  cigar   depth   strand_bias     rms_map_qual    in_hom_run      num_mapq_zero   num_alleles     num_reads_w_dels        haplotype_score qual_depth      allele_count    allele_bal      in_hm2  in_hm3  is_somatic      in_esp  aaf_esp_ea      aaf_esp_aa      aaf_esp_all     exome_chip      in_1kg  aaf_1kg_amr     aaf_1kg_asn     aaf_1kg_afr     aaf_1kg_eur     aaf_1kg_all     grc     gms_illumina    gms_solid       gms_iontorrent  in_cse  encode_tfbs     encode_dnaseI_cell_count        encode_dnaseI_cell_list encode_consensus_gm12878        encode_consensus_h1hesc encode_consensus_helas3 encode_consensus_hepg2  encode_consensus_huvec  encode_consensus_k562   gts     gt_types        gt_phases       gt_depths       gt_ref_depths   gt_alt_depths   gt_quals
    1       SMS173  1       chr1    100336360       100336361       60429   1       C       T       25701.56        None    snp     ts      1.0     1       rs2230306       None    None    None    None    None    None    None    None    None    None    None    None    chr1p21.2       None    0       0       1       None    2.24376e-65     2       6       4       0       0.583333333333  0.921158650238  -0.0285714285714        0.507246376812  0.274757        AGL     ENST00000361522 1       1       0       5       ctC/ctT L281    1515    protein_coding  synonymous_coding       LOW     None    None    None    None    None    None    None    1452    None    70.01   1       0       24      0.0     1.3604  19.85   14      None    None    None    None    1       0.304251        0.091728        0.232894        0       1       0.7     0.68    0.95    0.67    0.74    None    None    None    None    0       CEBPB_1 2       HCM;HCPEpiC     T       R       T       R       T       R       C|T,T|T,C|T,C||T,C|T,T|T,T|T,C|T,T|T,C|C,C|C,C|T 1,3,1,1,1,3,3,1,3,0,0,1 False,False,False,False,False,False,False,False,False,False,False,False 161,151,131,168,115,132,103,122,106,74,83,106   81,3,66,82,62,1,1,59,4,70,80,48 80,148,65,86,53,130,102,63,102,4,3,58   99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,80.05,99.0,99.0
    1       SMS173  1       chr1    100358102       100358103       60456   1       C       T       9734.77 None    snp     ts      1.0     1       rs3753494       None    None    None    None    None    None    None    None    None    None    None    GDE_C   chr1p21.2       None    0       0       1       None    2.26616e-55     8       3       1       0       0.208333333333  0.401650457515  0.242105263158  0.344202898551  0.243448        AGL     ENST00000361522 1       1       0       22      Cct|Tct P1050S  1515    protein_coding  non_syn_coding  MED     None    None    None    None    None    None    None    1476    None    70.03   0       0       24      0.0     1.8167  16.42   5       None    None    None    None    1       0.146163        0.126419        0.139474        1       1       0.12    0.02    0.14    0.15    0.11    None    None    None    None    0       None    None    None    T       R       T       T       R       T       C|T,C|C,C|C,C|C,C|C,T|T,C|T,C|T,C|C,C|C,C|C,C|C 1,0,0,0,0,3,1,1,0,0,0,0 False,False,False,False,False,False,False,False,False,False,False,False 213,122,152,169,114,143,119,118,106,69,55,96    108,119,152,166,113,7,59,64,104,67,53,92        105,3,0,3,1,136,60,54,2,2,2,4   99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0
    1       SMS173  2       chr1    15808197        15808198        14245   3       C       T       8880.77 None    snp     ts      1.0     1       rs7520335       None    None    None    None    None    None    None    None    None    None    None    None    chr1p36.21      None    0       1       0       None    None    7       5       0       0       0.208333333333  0.36197632685   -0.263157894737 0.344202898551  0.248348        CELA2B  ENST00000375909 1       1       0       4       Cgt/Tgt R68C    113     protein_coding  non_syn_coding  MED     None    None    None    None    None    None    None    1549    None    69.51   0   0       24      0.0     1.3894  12.7    5       None    None    None    None    0       None    None    None    0       1       0.22    0.53    0.19    0.25    0.31    None    None    None    None    0   None    None    None    R       R       T       R       T       R       C|T,C|C,C|T,C|C,C|T,C|C,C|T,C|T,C|C,C|C,C|C,C|C 1,0,1,0,1,0,1,1,0,0,0,0 False,False,False,False,False,False,False,False,False,False,False,False 214,134,199,233,86,172,83,117,91,55,61,104      125,131,111,231,50,171,28,62,91,53,61,104       89,3,88,2,36,0,55,55,0,2,0,0    99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,96.6,99.0,99.0
    1       SMS173  2       chr1    15808766        15808767        14249   2       G       A       3435.51 None    snp     ts      1.0     1       rs3820071       None    None    None    None    None    None    None    None    None    None    None    Trypsin chr1p36.21      None    0       1       0       None    6.64484e-08     7       5       0       0       0.208333333333  0.36197632685   -0.263157894737 0.344202898551  0.248209        CELA2B  ENST00000375910 1       1       0       4       Ggg/Agg G79R    269     protein_coding  non_syn_coding  MED     None    None    None    None    None    None    None    678     None    70.0    0       0       24      0.0     0.5304  11.08   5       None    None    None    None    1       0.245698        0.260781        0.250807        1       1       0.31    0.54    0.25    0.26    0.34    None    None    None    None    0       None    None    None    T       R       R       R       R       unknown G|A,G|G,G|A,G|G,G|A,G|G,G|A,G|A,G|G,G|G,G|G,G|G 1,0,1,0,1,0,1,1,0,0,0,0 False,False,False,False,False,False,False,False,False,False,False,False 86,53,101,106,50,58,35,38,46,25,34,46   55,51,57,104,32,55,19,23,45,25,33,46    31,2,44,2,18,3,16,15,1,0,1,0    99.0,56.03,99.0,99.0,99.0,76.93,99.0,99.0,91.93,69.16,59.8,99.0

This indicates that sample SMS173 has a candidate compound heterozygote in
AGL and CELA2B.

---------------------
``--columns``
---------------------

By default, this tool reports all columns in the ``variants`` table. One may
choose to report only a subset of the columns using the ``--columns`` option.  For
example, to report just the ``chrom, start, end, ref``, and ``alt`` columns, one
would use the following:

.. code-block:: bash

    $ gemini gemini comp_hets \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        my.db \
        | head -11

    family	sample	comp_het_id	gene	chrom	start	end	ref	alt	impact	impact_severity
    1	SMS173	1	AGL	chr1	100336360	100336361	C	T	synonymous_coding	LOW
    1	SMS173	1	AGL	chr1	100358102	100358103	C	T	non_syn_coding	MED
    1	SMS173	2	CELA2B	chr1	15808197	15808198	C	T	non_syn_coding	MED
    1	SMS173	2	CELA2B	chr1	15808766	15808767	G	A	non_syn_coding	MED
    1	SMS173	3	CELA2B	chr1	15808197	15808198	C	T	non_syn_coding	MED
    1	SMS173	3	CELA2B	chr1	15808871	15808872	G	A	non_syn_coding	MED
    1	SMS173	4	CELA2B	chr1	15808766	15808767	G	A	non_syn_coding	MED
    1	SMS173	4	CELA2B	chr1	15808871	15808872	G	A	non_syn_coding	MED
    1	SMS173	5	AJAP1	chr1	4772052	4772053	T	C	synonymous_coding	LOW
    1	SMS173	5	AJAP1	chr1	4834605	4834606	T	C	UTR_3_prime	LOW

.. note::

    The output will always start with the family ID, the sample name, and
    the compound heterozygote identification number.

--------------------
``--only-affected``
--------------------
By default, candidate compound heterozygous variants are reported for all
individuals in the database.  One can restrict the analysis to variants in
only individuals with an affected phenotype using the ``--only-affected`` option.

.. code-block:: bash

	$ gemini comp_hets --only-affected my.db

--------------------
``--families``
--------------------
By default, candidate compound heterozygous variants are reported for families
in the database.  One can restrict the analysis to variants in
specific familes with the ``--families`` option.  Families should be provided
as a comma-separated list

.. code-block:: bash

    $ gemini comp_hets --families 1 my.db
    $ gemini comp_hets --families 1,7 my.db


---------------------
``--ignore-phasing``
---------------------
If your genotypes aren't phased, we can't be certain that two heterozygotes
are on opposite alleles.  However, we can still identify pairs of heterozygotes
that are *candidates* for compound heterozygotes. Just use the
``--ignore-phasing`` option.


---------------------
``--filter``
---------------------

By default, this tool will report all variants regardless of their putative
functional impact.  In order to apply additional constraints on the variants
returned, one can use the ``--filter`` option.  Using SQL syntax, conditions
applied with the ``--filter option become WHERE clauses in the query issued to
the GEMINI database.  For example, if we wanted to restrict candidate variants
to solely those with a HIGH predicted functional consequence, we could use the
following:

.. code-block:: bash

    $ gemini gemini comp_hets \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --filter "impact_severity = 'HIGH'"
        my.db \
        | head -11

    family	sample	comp_het_id	gene	chrom	start	end	ref	alt	impact	impact_severity
    1	SMS173	1	TMCO4	chr1	20020993	20020994	C	CGT	frame_shift	HIGH
    1	SMS173	1	TMCO4	chr1	20020994	20020995	G	GTG	frame_shift	HIGH
    1	SMS173	2	HRNR	chr1	152185788	152185789	G	GCGACTAGG	frame_shift	HIGH
    1	SMS173	2	HRNR	chr1	152187906	152187907	T	TA	frame_shift	HIGH
    1	SMS173	3	FAM131C	chr1	16384996	16384997	G	GCA	frame_shift	HIGH
    1	SMS173	3	FAM131C	chr1	16384998	16384999	G	GCA	frame_shift	HIGH
    1	SMS173	4	CEP104	chr1	3753055	3753056	T	TTTTT	splice_donor	HIGH
    1	SMS173	4	CEP104	chr1	3753056	3753057	A	T	splice_donor	HIGH
    1	SMS173	5	AL355149.1	chr1	16862565	16862566	G	A	stop_gain	HIGH
    1	SMS173	5	AL355149.1	chr1	16863313	16863314	A	ACCCCTTTCTGCTG	frame_shift	HIGH




===========================================================================
``de_novo``: Identifying potential de novo mutations.
===========================================================================
.. note::

    1. This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:

    ``gemini load -v my.vcf -p my.ped my.db``


`Example PED file format for GEMINI`

.. code-block:: bash

	#Family_ID	Individual_ID	Paternal_ID	Maternal_ID	Sex	Phenotype	Ethnicity
	1	S173	S238	S239	1	2	caucasian
	1	S238	-9	-9	1	1	caucasian
	1	S239	-9	-9	2	1	caucasian
	2	S193	S230	S231	1	2	caucasian
	2	S230	-9	-9	1	1	caucasian
	2	S231	-9	-9	2	1	caucasian
	3	S242	S243	S244	1	2	caucasian
	3	S243	-9	-9	1	1	caucasian
	3	S244	-9	-9	2	1	caucasian
	4	S253	S254	S255	1	2	caucasianNEuropean
	4	S254	-9	-9	1	1	caucasianNEuropean
	4	S255	-9	-9	2	1	caucasianNEuropean


Assuming you have defined the familial relationships between samples when loading
your VCF into GEMINI, one can leverage a built-in tool for identifying de novo
(a.k.a spontaneous) mutations that arise in offspring.

---------------------
``default behavior``
---------------------

By default, the ``de novo`` tool will report, for each
family in the database, all columns in the variants table for mutations that
are not found in the parents yet are observed as heterozygotes in the offspring.
For example:

.. code-block:: bash

    $ gemini de_novo my.db

    family_id	family_members	family_genotypes   family_genotype_depths	chrom	start	end	variant_id	anno_id	ref	alt	qual	filter	type	sub_type	call_rate	in_dbsnp	rs_ids	in_omim	clinvar_sig	clinvar_disease_name	clinvar_dbsource	clinvar_dbsource_id	clinvar_origin	clinvar_dsdb	clinvar_dsdbid	clinvar_disease_acc	clinvar_in_locus_spec_db	clinvar_on_diag_assay	pfam_domain	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	gerp_bp_score	gerp_element_pval	num_hom_ref	num_het	num_hom_alt	num_unknown	aaf	hwe	inbreeding_coeff	pi	recomb_rate	gene	transcript	is_exonic	is_coding	is_lof	exon	codon_change	aa_change	aa_length	biotype	impact	impact_severity	polyphen_pred	polyphen_score	sift_pred	sift_score	anc_allele	rms_bq	cigar	depth	strand_bias	rms_map_qual	in_hom_run	num_mapq_zero	num_alleles	num_reads_w_dels	haplotype_score	qual_depth	allele_count	allele_bal	in_hm2	in_hm3	is_somatic	in_esp	aaf_esp_ea	aaf_esp_aa	aaf_esp_all	exome_chip	in_1kg	aaf_1kg_amr	aaf_1kg_asn	aaf_1kg_afr	aaf_1kg_eur	aaf_1kg_all	grc	gms_illumina	gms_solid	gms_iontorrent	in_cse	encode_tfbs	encode_dnaseI_cell_count	encode_dnaseI_cell_list	encode_consensus_gm12878	encode_consensus_h1hesc	encode_consensus_helas3	encode_consensus_hepg2	encode_consensus_huvec	encode_consensus_k562	gts	gt_types	gt_phases	gt_depths	gt_ref_depths	gt_alt_depths	gt_quals
    1	238(father; unknown),239(mother; unknown),173(child; affected)	AA/AA,AA/AA,AA/A	1,4,7	chr1	10067	10069	1	1	AA	A	113.21	None	indel	del	0.75	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Simple_repeat_Simple_repeat_(CCCTAA)n;trf;Satellite_telo_TAR1;trf;trf;trf;trf;trf	0	1	0	None	None	6	1	2	3	0.277777777778	0.0300651703342	0.723076923077	0.424836601307	2.981822	WASH7P	ENST00000423562	0	0	0	None	None	None	None	unprocessed_pseudogene	downstream	LOW	None	None	None	None	None	None	None	212	None	11.39	1	84	18	None	30.4532	1.55	5	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	91.7	47.1	94.7	0	None	None	None	CTCF	CTCF	unknown	unknown	unknown	CTCF	AA/A,./.,A/A,AA/AA,AA/AA,AA/AA,A/A,AA/AA,AA/AA,./.,AA/AA,./.	1,2,3,0,0,0,3,0,0,2,0,2	False,False,False,False,False,False,False,False,False,False,False,False	7,-1,2,4,1,4,2,2,1,-1,1,-1	33,-1,28,33,11,12,7,23,7,-1,12,-1	1,-1,2,0,0,0,2,0,0,-1,0,-1	26.74,-1.0,6.02,12.04,3.01,11.81,6.02,6.02,3.01,-1.0,3.01,-1.0
    4	254(father; unknown),255(mother; unknown),253(child; affected)	G/G,G/G,G/A	38,19,21	chr1	13109	13110	4	1	G	A	34.7	None	snp	ts	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	9	3	0	0	0.125	0.620690717057	-0.142857142857	0.228260869565	2.981822	WASH7	ENST00000423562	0	0	0	None	None	None	None	unprocessed_pseudogene	downstream	LOW	None	None	None	None	None	None	None	458	None	30.96	1	14	24	0.0	2.317	0.32	3	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	unknown	R	unknown	T	G/G,G/G,G/G,G/A,G/G,G/G,G/G,G/A,G/G,G/A,G/G,G/G	0,0,0,1,0,0,0,1,0,1,0,0	False,False,False,False,False,False,False,False,False,False,False,False	55,28,101,54,29,53,14,34,12,21,38,19	55,27,97,42,28,51,13,31,12,18,34,16	0,1,4,12,1,2,1,3,0,3,4,3	81.18,11.7,99.0,59.65,51.14,40.46,18.05,24.49,18.04,3.35,69.19,5.41
    1	238(father; unknown),239(mother; unknown),173(child; affected)	GTTG/GTTG,GTTG/GTTG,GTTG/G	21,59,41	chr1	14398	14402	13	1	GTTG	G	97.43	None	indel	del	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	9	3	0	0	0.125	0.620690717057	-0.142857142857	0.228260869565	2.981822	DDX11L1	ENST00000450305	0	0	0	None	None	None	None	transcribed_unprocessed_pseudogene	downstream	LOW	None	None	None	None	None	None	None	2045	None	15.9	0	4	24	None	145.8039	0.13	3	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	0.0	0.0	43.5	0	None	None	None	R	R	CTCF	R	R	T	GTTG/G,GTTG/G,GTTG/GTTG,GTTG/G,GTTG/GTTG,GTTG/GTTG,GTTG/GTTG,GTTG/GTTG,GTTG/GTTG,GTTG/GTTG,GTTG/GTTG,GTTG/GTTG	1,1,0,1,0,0,0,0,0,0,0,0	False,False,False,False,False,False,False,False,False,False,False,False	41,56,69,35,21,59,21,27,8,23,33,15	226,225,235,235,143,214,111,124,115,105,128,101	23,23,15,13,0,1,0,0,0,1,0,5	81.0,36.2,99.0,48.04,63.22,24.03,63.22,81.27,24.08,69.24,48.14,45.15
    1	238(father; unknown),239(mother; unknown),173(child; affected)	A/A,A/A,A/G	152,214,250	chr1	14541	14542	18	1	A	G	1369.37	None	snp	ts	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	4	8	0	0	0.333333333333	0.0832645169833	-0.5	0.463768115942	2.981822	DDX11L1	ENST00000456328	0	0	0	None	None	None	None	processed_transcript	downstream	LOW	None	None	None	None	None	None	None	2095	None	19.42	1	105	24	0.0	0.8894	1.01	8	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	CTCF	R	R	T	A/G,A/G,A/A,A/G,A/A,A/A,A/G,A/G,A/G,A/G,A/A,A/G	1,1,0,1,0,0,1,1,1,1,0,1	False,False,False,False,False,False,False,False,False,False,False,False	250,247,250,250,152,214,124,171,81,96,124,136	212,231,235,229,144,198,104,162,66,83,114,125	38,16,15,21,8,16,20,9,15,13,10,10	99.0,66.22,99.0,99.0,22.53,26.79,99.0,63.15,99.0,32.64,47.1,99.0    ...
    ...

.. note::

    The output will always start with the family ID, the family members, the
    observed genotypes, and the observed aligned sequencing depths
    for the family members.

---------------------
``-d [0]``
---------------------

Unfortunately, inherited variants can often appear to be de novo mutations simply because
insufficient sequence coverage was available for one of the parents to detect that the
parent(s) is also a heterozygote (and thus the variant was actually inherited, not
spontaneous).  One simple way to filter such artifacts is to enforce a minimum sequence
depth (default: 0) for each sample.  For example, if we require that at least 50 sequence
alignments were present for mom, dad and child, two of the above variants will be eliminated
as candidates:

.. code-block:: bash

    $ gemini de_novo -d 50 my.db

    family_id	family_members	family_genotypes   family_genotype_depths	chrom	start	end	variant_id	anno_id	ref	alt	qual	filter	type	sub_type	call_rate	in_dbsnp	rs_ids	in_omim	clinvar_sig	clinvar_disease_name	clinvar_dbsource	clinvar_dbsource_id	clinvar_origin	clinvar_dsdb	clinvar_dsdbid	clinvar_disease_acc	clinvar_in_locus_spec_db	clinvar_on_diag_assay	pfam_domain	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	gerp_bp_score	gerp_element_pval	num_hom_ref	num_het	num_hom_alt	num_unknown	aaf	hwe	inbreeding_coeff	pi	recomb_rate	gene	transcript	is_exonic	is_coding	is_lof	exon	codon_change	aa_change	aa_length	biotype	impact	impact_severity	polyphen_pred	polyphen_score	sift_pred	sift_score	anc_allele	rms_bq	cigar	depth	strand_bias	rms_map_qual	in_hom_run	num_mapq_zero	num_alleles	num_reads_w_dels	haplotype_score	qual_depth	allele_count	allele_bal	in_hm2	in_hm3	is_somatic	in_esp	aaf_esp_ea	aaf_esp_aa	aaf_esp_all	exome_chip	in_1kg	aaf_1kg_amr	aaf_1kg_asn	aaf_1kg_afr	aaf_1kg_eur	aaf_1kg_all	grc	gms_illumina	gms_solid	gms_iontorrent	in_cse	encode_tfbs	encode_dnaseI_cell_count	encode_dnaseI_cell_list	encode_consensus_gm12878	encode_consensus_h1hesc	encode_consensus_helas3	encode_consensus_hepg2	encode_consensus_huvec	encode_consensus_k562	gts	gt_types	gt_phases	gt_depths	gt_ref_depths	gt_alt_depths	gt_quals
    1	238(father; unknown),239(mother; unknown),173(child; affected)	A/A,A/A,A/G	152,214,250	chr1	14541	14542	18	1	A	G	1369.37	None	snp	ts	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	4	8	0	0	0.333333333333	0.0832645169833	-0.5	0.463768115942	2.981822	DDX11L1	ENST00000456328	0	0	0	None	None	None	None	processed_transcript	downstream	LOW	None	None	None	None	None	None	None	2095	None	19.42	1	105	24	0.0	0.8894	1.01	8	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	CTCF	R	R	T	A/G,A/G,A/A,A/G,A/A,A/A,A/G,A/G,A/G,A/G,A/A,A/G	1,1,0,1,0,0,1,1,1,1,0,1	False,False,False,False,False,False,False,False,False,False,False,False	250,247,250,250,152,214,124,171,81,96,124,136	212,231,235,229,144,198,104,162,66,83,114,125	38,16,15,21,8,16,20,9,15,13,10,10	99.0,66.22,99.0,99.0,22.53,26.79,99.0,63.15,99.0,32.64,47.1,99.0
    1	238(father; unknown),239(mother; unknown),173(child; affected)	A/A,A/A,A/G	189,250,250	chr1	14573	14574	19	1	A	G	723.72	None	snp	ts	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	6	6	0	0	0.25	0.248213079014	-0.333333333333	0.391304347826	2.981822	DDX11L1	ENST00000456328	0	0	0	None	None	None	None	processed_transcript	downstream	LOW	None	None	None	None	None	None	None	2233	None	20.21	0	73	24	0.0	1.1058	0.63	6	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	CTCF	R	R	T	A/G,A/G,A/A,A/G,A/A,A/A,A/G,A/G,A/G,A/A,A/A,A/A	1,1,0,1,0,0,1,1,1,0,0,0	False,False,False,False,False,False,False,False,False,False,False,False	250,248,250,241,189,250,130,189,92,107,146,141	218,232,237,221,181,232,115,177,76,97,136,134	32,14,13,20,8,17,15,12,16,10,10,7	99.0,31.97,99.0,99.0,96.41,99.0,64.51,35.62,99.0,26.4,65.9,0.76
    1	238(father; unknown),239(mother; unknown),173(child; affected)	G/G,G/G,G/A	197,247,250	chr1	14589	14590	20	1	G	A	178.22	None	snp	ts	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	8	4	0	0	0.166666666667	0.488422316764	-0.2	0.289855072464	2.981822	DDX11L1	ENST00000456328	0	0	0	None	None	None	None	processed_transcript	downstream	LOW	None	None	None	None	None	None	None	2234	None	21.45	0	37	24	0.0	0.9191	0.25	4	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	CTCF	R	R	T	G/A,G/G,G/G,G/A,G/G,G/G,G/A,G/G,G/A,G/G,G/G,G/G	1,0,0,1,0,0,1,0,1,0,0,0	False,False,False,False,False,False,False,False,False,False,False,False	250,238,250,233,197,247,134,192,97,109,149,137	227,228,239,213,186,227,124,181,84,105,144,128	23,10,11,20,11,20,10,11,13,4,5,9	99.0,99.0,99.0,25.64,99.0,99.0,31.54,19.87,54.49,97.64,99.0,42.52
    1	238(father; unknown),239(mother; unknown),173(child; affected)	T/T,T/T,T/A	195,250,249	chr1	14598	14599	21	1	T	A	44.09	None	snp	tv	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	10	2	0	0	0.0833333333333	0.752823664836	-0.0909090909091	0.159420289855	2.981822	DDX11L1	ENST00000456328	0	0	0	None	None	None	None	processed_transcript	downstream	LOW	None	None	None	None	None	None	None	2245	None	22.1	0	18	24	0.0	1.1988	0.13	2	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	CTCF	R	R	T	T/A,T/T,T/T,T/T,T/T,T/T,T/T,T/T,T/A,T/T,T/T,T/T	1,0,0,0,0,0,0,0,1,0,0,0	False,False,False,False,False,False,False,False,False,False,False,False	249,237,250,242,195,250,138,209,91,102,148,133	226,229,240,223,187,231,129,198,76,94,140,118	23,8,10,19,8,19,9,11,15,8,8,14	65.38,99.0,99.0,92.74,99.0,99.0,23.58,84.54,30.04,99.0,99.0,45.7
    ...


---------------------
``--columns``
---------------------

By default, this tool reports all columns in the ``variants`` table. One may
choose to report only a subset of the columns using the ``--columns`` option.  For
example, to report just the ``chrom, start, end, ref``, and ``alt`` columns, one
would use the following:

.. code-block:: bash

    $ gemini de_novo -d 50 --columns "chrom, start, end, ref, alt" my.db

    family_id	family_members	family_genotypes   family_genotype_depths	chrom	start	end	ref	alt
    1	238(father; unknown),239(mother; unknown),173(child; affected)	A/A,A/A,A/G	152,214,250	chr1	14541	14542	A	G
    1	238(father; unknown),239(mother; unknown),173(child; affected)	A/A,A/A,A/G	189,250,250	chr1	14573	14574	A	G
    1	238(father; unknown),239(mother; unknown),173(child; affected)	G/G,G/G,G/A	197,247,250	chr1	14589	14590	G	A
    1	238(father; unknown),239(mother; unknown),173(child; affected)	T/T,T/T,T/A	195,250,249	chr1	14598	14599	T	A
    ...

.. note::

    The output will always start with the family ID, the family members, the
    observed genotypes, and the observed aligned sequencing depths
    for the family members.


---------------------
``--filter``
---------------------

By default, this tool will report all variants regardless of their putative
functional impact.  In order to apply additional constraints on the variants
returned, one can use the ``--filter`` option.  Using SQL syntax, conditions
applied with the ``--filter option become WHERE clauses in the query issued to
the GEMINI database.  For example, if we wanted to restrict candidate variants
to solely those with a HIGH predicted functional consequence, we could use the
following:

.. code-block:: bash

    $ gemini de_novo -d 50 \
          --columns "chrom, start, end, ref, alt" \
          --filter "impact_severity = 'HIGH'" \
          my.db

    family_id	family_members	family_genotypes   family_genotype_depths	chrom	start	end	ref	alt
    3	243(father; unknown),244(mother; unknown),242(child; affected)	C/C,C/C,C/A	249,243,250	chr1	17729	17730	C	A
    4	254(father; unknown),255(mother; unknown),253(child; affected)	A/A,A/A,A/G	86,146,83	chr1	168097	16809	A	G
    4	254(father; unknown),255(mother; unknown),253(child; affected)	G/G,G/G,G/T	107,182,72	chr1	12854400	12854401	G	T
    3	243(father; unknown),244(mother; unknown),242(child; affected)	A/A,A/A,A/ATGGTGTTG	211,208,208	chr1	12855995	12855996	A	ATGGTGTTG
    ...

-------------------------
``--min-kindreds [None]``
-------------------------
By default, the ``de_novo`` tool will report every candidate mutation variant
that impacts at least one of the families in the database.  However, one
can restrict the reported genes to those where de novo variants
were observed in more than one family (thus further substantiating the potential role of the gene in the etiology of the phenotype).

For example, the following command would further restrict candidate genes to those genes with a de novo variant in at least two families:

.. code-block:: bash

    $ gemini de_novo -d 50 \
          --columns "chrom, start, end, ref, alt" \
          --filter "impact_severity = 'HIGH'" \
          --min-kindreds 2 \
          my.db


--------------------
``--only-affected``
--------------------
By default, candidate de novo mutations are reported for all
individuals in the database.  One can restrict the analysis to variants in
only individuals with an affected phenotype using the ``--only-affected`` option.

.. code-block:: bash

    $ gemini de_novo --only-affected my.db


--------------------
``--families``
--------------------
By default, candidate de novo variants are reported for families
in the database.  One can restrict the analysis to variants in
specific familes with the ``--families`` option.  Families should be provided
as a comma-separated list

.. code-block:: bash

    $ gemini de_novo --families 1 my.db
    $ gemini de_novo --families 1,7 my.db


============================================================================
``autosomal_recessive``: Find variants meeting an autosomal recessive model.
============================================================================
.. warning::

    By default, this tool requires that you identify familial relationships
    via a PED file when loading your VCF into GEMINI.  For example:

    ``gemini load -v my.vcf -p my.ped my.db``

    However, in the absence of established parent/child relationships in the PED
    file, GEMINI will issue a WARNING, yet will attempt to identify autosomal
    recessive candidates for all samples marked as "affected".

---------------------
``default behavior``
---------------------

Assuming you have defined the familial relationships between samples when
loading your VCF into GEMINI, one can leverage a built-in tool for
identifying variants that meet an autosomal recessive inheritance pattern.
The reported variants will be restricted to those variants having the
potential to impact the function of affecting protein coding transcripts.

For the following examples, let's assume we have a PED file for 3 different
families as follows (the kids are affected in each family, but the parents
are not):

.. code-block:: bash

    $ cat families.ped
    1	1_dad	0	0	-1	1
    1	1_mom	0	0	-1	1
    1	1_kid	1_dad	1_mom	-1	2
    2	2_dad	0	0	-1	1
    2	2_mom	0	0	-1	1
    2	2_kid	2_dad	2_mom	-1	2
    3	3_dad	0	0	-1	1
    3	3_mom	0	0	-1	1
    3	3_kid	3_dad	3_mom	-1	2

.. code-block:: bash

    $ gemini autosomal_recessive my.db
    family_id   family_members  family_genotypes    family_genotype_depths  chrom   start   end variant_id  anno_id ref alt qual    filter  type    sub_type    call_rate   in_dbsnp    rs_ids  in_omim clinvar_sig clinvar_disease_name    clinvar_dbsource    clinvar_dbsource_id clinvar_origin  clinvar_dsdb    clinvar_dsdbid  clinvar_disease_acc clinvar_in_locus_spec_db    clinvar_on_diag_assay   pfam_domain cyto_band   rmsk    in_cpg_island   in_segdup   is_conserved    gerp_bp_score   gerp_element_pval   num_hom_ref num_het num_hom_alt num_unknown aaf hwe inbreeding_coeff    pi  recomb_rate gene    transcript  is_exonic   is_coding   is_lof  exon    codon_change    aa_change   aa_length   biotype impact  impact_severity polyphen_pred   polyphen_score  sift_pred   sift_score  anc_allele  rms_bq  cigar   depth   strand_bias rms_map_qual    in_hom_run  num_mapq_zero   num_alleles num_reads_w_dels    haplotype_score qual_depth  allele_count    allele_bal  in_hm2  in_hm3  is_somatic  in_esp  aaf_esp_ea  aaf_esp_aa  aaf_esp_all exome_chip  in_1kg  aaf_1kg_amr aaf_1kg_asn aaf_1kg_afr aaf_1kg_eur aaf_1kg_all grc gms_illumina    gms_solid   gms_iontorrent  in_cse  encode_tfbs encode_dnaseI_cell_count    encode_dnaseI_cell_list encode_consensus_gm12878    encode_consensus_h1hesc encode_consensus_helas3 encode_consensus_hepg2  encode_consensus_huvec  encode_consensus_k562   gts gt_types    gt_phases   gt_depths   gt_ref_depths   gt_alt_depths   gt_quals
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  C/T,C/T,T/T 39,29,24    chr10   48004991    48004992    3   1   C   T   1047.87 None    snp ts  1.0 0   None    None    None    None    None    None    None    None    None    None    None    None    None    chr10q11.22 None    0   1   0   None    None    0   8   1   0   0.555555555556  0.0163950703837 -0.8    0.522875816993  1.718591    ASAH2C  ENST00000420079 1   1   0   exon_10_48003968_48004056   tGt/tAt C540Y   610 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    165 None    20.94   0   0   8   0.0 4.383   9.53    4   None    None    None    None    0   None    None    Non 0   0   None    None    None    None    None    grc_fix None    None    None    0   None    None    None    R   R   R   R   R   R   C/T,C/T,C/T,C/T,C/T,T/T,C/T,C/T,C/T 1,1,1,1,1,3,1,1,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  C/T,C/T,T/T 39,29,24    chr10   48003991    48003992    2   1   C   T   1047.87 None    snp ts  1.0 1   rs142685947 None    None    None    None    None    None    None    None    None    None    None    None    chr10q11.22 None    0   1   1   None    3.10871e-42 0   8   1   0   0.555555555556  0.0163950703837 -0.8    0.522875816993  1.718591    ASAH2C  ENST00000420079 1   1   0   exon_10_48003968_48004056   tGt/tAt C540Y   610 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    165 None    20.94   0   0   8   0.0 4.383   9.53    4   None    None    None    None    0   Non None    None    0   0   None    None    None    None    None    grc_fix 73.3    40.3    92.8    0   None    None    None    R   R   R   R   R   R   C/T,C/T,T/T,C/T,C/T,C/T,C/T,C/T,C/T 1,1,3,1,1,1,1,1,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    3   3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)  T/C,T/C,C/C 39,29,24    chr10   135369531   135369532   5   6   T   C   122.62  None    snp ts  1.0 1   rs3747881   None    None    None    None    None    None    None    None    None    None    None    None    chr10q26.3  None    0   0   1   None    3.86096e-59 0   8   1   0   0.555555555556  0.0163950703837 -0.8    0.522875816993  0.022013    SYCE1   ENST00000368517 1   1   0   exon_10_135369485_135369551 aAg/aGg K147R   282 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    239 None    36.02   2   0   8   0.0 5.7141  2.31    2   None    None    None    None    1   0.093837    0.163867    0.117561    1   0   None    None    None    None    None    None    None    None    None    0   None    None    None    R   R   R   R   R   R   T/C,T/C,T/C,T/C,T/C,T/C,T/C,T/C,C/C 1,1,1,1,1,1,1,1,3   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  T/C,T/C,C/C 39,29,24    chr10   1142207 1142208 1   4   T   C   3404.3  None    snp ts  1.0 1   rs10794716  None    None    None    None    None    None    None    None    None    None    None    None    chr10p15.3  None    0   0   0   None    None    0   7   2   0   0.611111111111  0.0562503650686 -0.636363636364 0.503267973856  0.200924    WDR37   ENST00000381329 1   1   1   exon_10_1142110_1142566 Tga/Cga *250R   249 protein_coding  stop_loss   HIG None    None    None    None    None    None    None    122 None    36.0    0   0   8   0.0 2.6747  27.9    8   None    None    None    None    1   0.000465    0.024966    0.008765    0   1   1   1   0.98    1   0.99    None    None    None    None    0   None    2   Osteobl;Progfib T   T   T   T   T   T   T/C,T/C,C/C,T/C,T/C,C/C,T/C,T/C,T/C 1,1,3,1,1,3,1,1,1   False,False,False,False,False,False,False,False,False   39,29,24,59,49,64,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  T/C,T/C,C/C 59,49,64    chr10   1142207 1142208 1   4   T   C   3404.3  None    snp ts  1.0 1   rs10794716  None    None    None    None    None    None    None    None    None    None    None    None    chr10p15.3  None    0   0   0   None    None    0   7   2   0   0.611111111111  0.0562503650686 -0.636363636364 0.503267973856  0.200924    WDR37   ENST00000381329 1   1   1   exon_10_1142110_1142566 Tga/Cga *250R   249 protein_coding  stop_loss   HIG None    None    None    None    None    None    None    122 None    36.0    0   0   8   0.0 2.6747  27.9    8   None    None    None    None    1   0.000465    0.024966    0.008765    0   1   1   1   0.98    1   0.99    None    None    None    None    0   None    2   Osteobl;Progfib T   T   T   T   T   T   T/C,T/C,C/C,T/C,T/C,C/C,T/C,T/C,T/C 1,1,3,1,1,3,1,1,1   False,False,False,False,False,False,False,False,False   39,29,24,59,49,64,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    ...

.. note::

    The output will always start with the family ID, the family members, the
    observed genotypes, and the observed aligned sequencing depths
    for the family members.

---------------------
``--columns``
---------------------

By default, this tool reports all columns in the ``variants`` table. One may
choose to report only a subset of the columns using the ``--columns`` option.  For
example, to report just the ``gene, chrom, start, end, ref, alt, impact``, and ``impact_severity`` columns, one
would use the following:

.. code-block:: bash

    $ gemini autosomal_recessive \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        my.db

    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  C/T,C/T,T/T 39,29,24    ASAH2C  chr10   48004991    48004992    C   T   non_syn_coding  MED
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  C/T,C/T,T/T 39,29,24    ASAH2C  chr10   48003991    48003992    C   T   non_syn_coding  MED
    3   3_dad(father; unaffected),3_mom(mother; unaffected),3_kid(child; affected)  T/C,T/C,C/C 39,29,24    SYCE1   chr10   135369531   135369532   T   C   non_syn_coding  MED
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  T/C,T/C,C/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  T/C,T/C,C/C 59,49,64    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH


----------------------
``--min-kindreds [1]``
----------------------
By default, the ``autosomal_recessive`` tool will report every gene variant
that impacts at least one of the families in the database.  However, one
can restrict the reported genes to those where autosomal recessive variants
were observed in more than one family (thus further substantiating the potential role of the gene in the etiology of the phenotype).

For example, to restricted the report to genes with variants (doesn't have
to be the _same_ variant) observed in at least two kindreds, use the following:


.. code-block:: bash

    $ gemini autosomal_recessive \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --min-kindreds 2 \
        my.db
    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  C/T,C/T,T/T 39,29,24    ASAH2C  chr10   48004991    48004992    C   T   non_syn_coding  MED
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  C/T,C/T,T/T 39,29,24    ASAH2C  chr10   48003991    48003992    C   T   non_syn_coding  MED
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  T/C,T/C,C/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  T/C,T/C,C/C 59,49,64    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH

--------------------
``--families``
--------------------
By default, candidate autosomal recessive variants are reported for families
in the database.  One can restrict the analysis to variants in
specific familes with the ``--families`` option.  Families should be provided
as a comma-separated list

.. code-block:: bash

    $ gemini autosomal_recessive --families 1 my.db
    $ gemini autosomal_recessive --families 1,7 my.db


---------------------
``--filter``
---------------------

By default, this tool will report all variants regardless of their putative
functional impact.  In order to apply additional constraints on the variants
returned, one can use the ``--filter`` option.  Using SQL syntax, conditions
applied with the ``--filter option become WHERE clauses in the query issued to
the GEMINI database.  For example, if we wanted to restrict candidate variants
to solely those with a HIGH predicted functional consequence, we could use the
following:

.. code-block:: bash

    $ gemini autosomal_recessive \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --min-kindreds 2 \
        --filter "impact_severity = 'HIGH'" \
        my.db

    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity
    1   1_dad(father; unaffected),1_mom(mother; unaffected),1_kid(child; affected)  T/C,T/C,C/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    2   2_dad(father; unaffected),2_mom(mother; unaffected),2_kid(child; affected)  T/C,T/C,C/C 59,49,64    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH

---------------------
``-d [0]``
---------------------

In order to eliminate less confident genotypes, it is possible to enforce a minimum sequence
depth (default: 0) for each sample:

.. code-block:: bash

    $ gemini autosomal_dominant \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --filter "impact_severity = 'HIGH'" \
        --min-kindreds 1 \
        -d 40 \
        my.db

    family_id   family_members  family_genotypes    gene    chrom   start   end ref alt impact  impact_severity
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    T/T,T/C,T/C WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   T/C,T/T,T/C WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH




===========================================================================
``autosomal_dominant``: Find variants meeting an autosomal dominant model.
===========================================================================

.. warning::

    1. By default, this tool requires that you identify familial relationships
    via a PED file when loading your VCF into GEMINI.  For example:

    ``gemini load -v my.vcf -p my.ped my.db``

    2. However, if neither parent is known to be affected, this tool will report any
       variant where one and only of the parents is heterozygous and the affected
       child is also heterozygous.  If one and only one of the parents is affected,
       the tool will report variants where both the affected child and the affected
       parent are heterozygous.  If both parents are known to be affected, the
       tool will report nothing for that family.  If parents are unknown, the tool
       will report variants where an affected individual is heterozygous and
       all unaffected individuals are homozygous for the reference allele.

---------------------
``default behavior``
---------------------

Assuming you have defined the familial relationships between samples when loading
your VCF into GEMINI, one can leverage a built-in tool for identifying variants
that meet an autosomal dominant inheritance pattern. The reported variants
will be restricted to those variants having the potential to impact the
function of affecting protein coding transcripts.

For the following examples, let's assume we have a PED file for 3 different
families as follows (the kids are affected in each family, but the parents
are not):

.. code-block:: bash

    $ cat families.ped
    1	1_dad	0	0	-1	1
    1	1_mom	0	0	-1	1
    1	1_kid	1_dad	1_mom	-1	2
    2	2_dad	0	0	-1	1
    2	2_mom	0	0	-1	2
    2	2_kid	2_dad	2_mom	-1	2
    3	3_dad	0	0	-1	2
    3	3_mom	0	0	-1	-9
    3	3_kid	3_dad	3_mom	-1	2


.. code-block:: bash

    $ gemini autosomal_dominant my.db | head

    family_id   family_members  family_genotypes    family_genotype_depths  chrom   start   end variant_id  anno_id ref alt qual    filter  type    sub_type    call_rate   in_dbsnp    rs_ids  in_omim clinvar_sig clinvar_disease_name    clinvar_dbsource    clinvar_dbsource_id clinvar_origin  clinvar_dsdb    clinvar_dsdbid  clinvar_disease_acc clinvar_in_locus_spec_db    clinvar_on_diag_assay   pfam_domain cyto_band   rmsk    in_cpg_island   in_segdup   is_conserved    gerp_bp_score   gerp_element_pval   num_hom_ref num_het num_hom_alt num_unknown aaf hwe inbreeding_coeff    pi  recomb_rate gene    transcript  is_exonic   is_coding   is_lof  exon    codon_change    aa_change   aa_length   biotype impact  impact_severity polyphen_pred   polyphen_score  sift_pred   sift_score  anc_allele  rms_bq  cigar   depth   strand_bias rms_map_qual    in_hom_run  num_mapq_zero   num_alleles num_reads_w_dels    haplotype_score qual_depth  allele_count    allele_bal  in_hm2  in_hm3  is_somatic  in_esp  aaf_esp_ea  aaf_esp_aa  aaf_esp_all exome_chip  in_1kg  aaf_1kg_amr aaf_1kg_asn aaf_1kg_afr aaf_1kg_eur aaf_1kg_all grc gms_illumina    gms_solid   gms_iontorrent  in_cse  encode_tfbs encode_dnaseI_cell_count    encode_dnaseI_cell_list encode_consensus_gm12878    encode_consensus_h1hesc encode_consensus_helas3 encode_consensus_hepg2  encode_consensus_huvec  encode_consensus_k562   gts gt_types    gt_phases   gt_depths   gt_ref_depths   gt_alt_depths   gt_quals
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   C/T,C/C,C/T 39,29,24    chr10   48003991    48003992    3   1   C   T   1047.87 None    snp ts  1.0 1   rs142685947 None    None    None    None    None    None    None    None    None    None    None    None    chr10q11.22 None    0   1   1   None    3.10871e-42 4   5   0   0   0.277777777778  0.248563248239  -0.384615384615 0.424836601307  1.718591    ASAH2C  ENST00000420079 1   1   0   exon_10_48003968_48004056   tGt/tAt C540Y   610 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    165 None    20.94   0   0   8   0.0 4.383   9.53    4   None    None    None    None    0   Non None    None    0   0   None    None    None    None    None    grc_fix 73.3    40.3    92.8    0   None    None    None    R   R   R   R   R   R   C/C,C/C,C/T,C/C,C/T,C/T,C/T,C/C,C/T 0,0,1,0,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   C/T,C/C,C/T 39,29,24    chr10   48004991    48004992    4   1   C   T   1047.87 None    snp ts  1.0 0   None    None    None    None    None    None    None    None    None    None    None    None    None    chr10q11.22 None    0   1   0   None    None    4   5   0   0   0.277777777778  0.248563248239  -0.384615384615 0.424836601307  1.718591    ASAH2C  ENST00000420079 1   1   0   exon_10_48003968_48004056   tGt/tAt C540Y   610 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    165 None    20.94   0   0   8   0.0 4.383   9.53    4   None    None    None    None    0   None    None    Non 0   0   None    None    None    None    None    grc_fix None    None    None    0   None    None    None    R   R   R   R   R   R   C/C,C/C,C/T,C/C,C/T,C/T,C/T,C/C,C/T 0,0,1,0,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    C/C,C/T,C/T 39,29,24    chr10   48003991    48003992    3   1   C   T   1047.87 None    snp ts  1.0 1   rs142685947 None    None    None    None    None    None    None    None    None    None    None    None    chr10q11.22 None    0   1   1   None    3.10871e-42 4   5   0   0   0.277777777778  0.248563248239  -0.384615384615 0.424836601307  1.718591    ASAH2C  ENST00000420079 1   1   0   exon_10_48003968_48004056   tGt/tAt C540Y   610 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    165 None    20.94   0   0   8   0.0 4.383   9.53    4   None    None    None    None    0   None    None    None    0   0   None    None    None    None    None    grc_fix 73.3    40.3    92.8    0   None    None    None    R   R   R   R   R   R   C/C,C/C,C/T,C/C,C/T,C/T,C/T,C/C,C/T 0,0,1,0,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    C/C,C/T,C/T 39,29,24    chr10   48004991    48004992    4   1   C   T   1047.87 None    snp ts  1.0 0   None    None    None    None    None    None    None    None    None    None    None    None    None    chr10q11.22 None    0   1   0   None    None    4   5   0   0   0.277777777778  0.248563248239  -0.384615384615 0.424836601307  1.718591    ASAH2C  ENST00000420079 1   1   0   exon_10_48003968_48004056   tGt/tAt C540Y   610 protein_coding  non_syn_coding  MED None    None    None    None    None    None    None    165 None    20.94   0   0   8   0.0 4.383   9.53    4   None    None    None    None    0   None    Non None    0   0   None    None    None    None    None    grc_fix None    None    None    0   None    None    None    R   R   R   R   R   R   C/C,C/C,C/T,C/C,C/T,C/T,C/T,C/C,C/T 0,0,1,0,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   G/A,G/G,G/A 39,29,24    chr10   135336655   135336656   5   1   G   A   38.34   None    snp ts  1.0 1   rs6537611   None    None    None    None    None    None    None    None    None    None    None    None    chr10q26.3  None    0   0   0   None    None    1   8   0   0   0.444444444444  0.0163950703837 -0.8    0.522875816993  0.43264 SPRN    ENST00000541506 0   0   0   None    None    None    151 protein_coding  intron  LOW None    None    None    Non None    None    None    2   None    37.0    4   0   4   0.0 0.0 19.17   4   None    None    None    None    0   None    None    None    0   0   None    None    None    Non None    None    None    None    None    0   None    None    None    R   R   R   R   unknown R   G/A,G/A,G/A,G/A,G/A,G/A,G/A,G/G,G/A 1,1,1,1,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    T/T,T/C,T/C 39,29,24    chr10   1142207 1142208 1   4   T   C   3404.3  None    snp ts  1.0 1   rs10794716  None    None    None    None    None    None    None    None    None    None    None    None    chr10p15.3  None    0   0   0   None    None    4   5   0   0   0.277777777778  0.248563248239  -0.384615384615 0.424836601307  0.200924    WDR37   ENST00000381329 1   1   1   exon_10_1142110_1142566 Tga/Cga *250R   249 protein_coding  stop_loss   HIG None    None    None    None    None    None    None    122 None    36.0    0   0   8   0.0 2.6747  27.9    8   None    None    None    None    1   0.000465    0.024966    0.008765    0   1   1   1   0.98    1   0.99    None    None    None    None    0   None    2   Osteobl;Progfib T   T   T   T   T   T   T/T,T/T,T/C,T/T,T/C,T/C,T/C,T/T,T/C 0,0,1,0,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   T/C,T/T,T/C 39,29,24    chr10   1142207 1142208 1   4   T   C   3404.3  None    snp ts  1.0 1   rs10794716  None    None    None    None    None    None    None    None    None    None    None    None    chr10p15.3  None    0   0   0   None    None    4   5   0   0   0.277777777778  0.248563248239  -0.384615384615 0.424836601307  0.200924    WDR37   ENST00000381329 1   1   1   exon_10_1142110_1142566 Tga/Cga *250R   249 protein_coding  stop_loss   HIG None    None    None    None    None    None    None    122 None    36.0    0   0   8   0.0 2.6747  27.9    8   None    None    None    None    1   0.000465    0.024966    0.008765    0   1   1   1   0.98    1   0.99    None    None    None    None    0   None    2   Osteobl;Progfib T   T   T   T   T   T   T/T,T/T,T/C,T/T,T/C,T/C,T/C,T/T,T/C 0,0,1,0,1,1,1,0,1   False,False,False,False,False,False,False,False,False   39,29,24,39,29,24,39,29,24  1,0,0,1,0,0,1,0,0   37,29,24,37,29,24,37,29,24  87.16,78.2,66.14,87.16,78.2,66.14,87.16,78.2,66.14


---------------------
``--columns``
---------------------

By default, this tool reports all columns in the ``variants`` table. One may
choose to report only a subset of the columns using the ``--columns`` option.  For
example, to report just the ``gene, chrom, start, end, ref, alt, impact``, and ``impact_severity`` columns, one
would use the following:

.. code-block:: bash

    $ gemini autosomal_dominant \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        my.db

    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   C/T,C/C,C/T 39,29,24    ASAH2C  chr10   48003991    48003992    C   T   non_syn_coding  MED
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   C/T,C/C,C/T 39,29,24    ASAH2C  chr10   48004991    48004992    C   T   non_syn_coding  MED
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    C/C,C/T,C/T 39,29,24    ASAH2C  chr10   48003991    48003992    C   T   non_syn_coding  MED
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    C/C,C/T,C/T 39,29,24    ASAH2C  chr10   48004991    48004992    C   T   non_syn_coding  MED
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   G/A,G/G,G/A 39,29,24    SPRN    chr10   135336655   135336656   G   A   intron  LOW
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    T/T,T/C,T/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   T/C,T/T,T/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH

.. note::

    The output will always start with the family ID, the family members, and the
    observed genotypes for the family members.


----------------------
``--min-kindreds [1]``
----------------------
By default, the ``autosomal_dominant`` tool will report every gene variant
that impacts at least one of the families in the database.  However, one
can restrict the reported genes to those where autosomal dominant variants
were observed in more than one family (thus further substantiating the potential role of the gene in the etiology of the phenotype).

For example, to restricted the report to genes with variants (doesn't have
to be the _same_ variant) observed in at least two kindreds, use the following:


.. code-block:: bash

    $ gemini autosomal_dominant \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --min-kindreds 2 \
        my.db

    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   C/T,C/C,C/T 39,29,24    ASAH2C  chr10   48003991    48003992    C   T   non_syn_coding  MED
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   C/T,C/C,C/T 39,29,24    ASAH2C  chr10   48004991    48004992    C   T   non_syn_coding  MED
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    C/C,C/T,C/T 39,29,24    ASAH2C  chr10   48003991    48003992    C   T   non_syn_coding  MED
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    C/C,C/T,C/T 39,29,24    ASAH2C  chr10   48004991    48004992    C   T   non_syn_coding  MED
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    T/T,T/C,T/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   T/C,T/T,T/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH


--------------------
``--families``
--------------------
By default, candidate autosomal dominant variants are reported for families
in the database.  One can restrict the analysis to variants in
specific familes with the ``--families`` option.  Families should be provided
as a comma-separated list

.. code-block:: bash

    $ gemini autosomal_dominant --families 1 my.db
    $ gemini autosomal_dominant --families 1,7 my.db


---------------------
``--filter``
---------------------

By default, this tool will report all variants regardless of their putative
functional impact.  In order to apply additional constraints on the variants
returned, one can use the ``--filter`` option.  Using SQL syntax, conditions
applied with the ``--filter option become WHERE clauses in the query issued to
the GEMINI database.  For example, if we wanted to restrict candidate variants
to solely those with a HIGH predicted functional consequence, we could use the
following:

.. code-block:: bash

    $ gemini autosomal_dominant \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --filter "impact_severity = 'HIGH'" \
        --min-kindreds 2 \
        my.db

    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity
    2   2_dad(father; unaffected),2_mom(mother; affected),2_kid(child; affected)    T/T,T/C,T/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH
    3   3_dad(father; affected),3_mom(mother; unknown),3_kid(child; affected)   T/C,T/T,T/C 39,29,24    WDR37   chr10   1142207 1142208 T   C   stop_loss   HIGH


---------------------
``-d [0]``
---------------------

In order to eliminate less confident genotypes, it is possible to enforce a minimum sequence
depth (default: 0) for each sample (in this case, no variants would meet this criteria):

.. code-block:: bash

    $ gemini autosomal_dominant \
        --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
        --filter "impact_severity = 'HIGH'" \
        --min-kindreds 1 \
        -d 40 \
        my.db

    family_id   family_members  family_genotypes    family_genotype_depths  gene    chrom   start   end ref alt impact  impact_severity



===========================================================================
``pathways``: Map genes and variants to KEGG pathways.
===========================================================================
Mapping genes to biological pathways is useful in understanding the
function/role played by a gene. Likewise, genes involved in common pathways
is helpful in understanding heterogeneous diseases. We have integrated
the KEGG pathway mapping for gene variants, to explain/annotate variation.
This requires your VCF be annotated with either snpEff/VEP.

Examples:

.. code-block:: bash

	$ gemini pathways -v 68 example.db
	chrom	start	end	ref	alt	impact	sample	genotype	gene	transcript	pathway
	chr10	52004314	52004315	T	C	intron	M128215	C/C	ASAH2	ENST00000395526	hsa00600:Sphingolipid_metabolism,hsa01100:Metabolic_pathways
	chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000531469	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer
	chr16	72057434	72057435	C	T	non_syn_coding	M10475	C/T	DHODH	ENST00000219240	hsa01100:Metabolic_pathways,hsa00240:Pyrimidine_metabolism


Here, -v specifies the version of the Ensembl genes used to build the KEGG
pathway map. Hence, use versions that match the VEP/snpEff versions of the
annotated vcf for correctness. For e.g VEP v2.6 and snpEff v3.1 use Ensembl
68 version of the genomes.

We currently support versions 66 through 71 of the Ensembl genes


---------------
``--lof``
---------------
By default, all gene variants that map to pathways are reported.  However,
one may want to restrict the analysis to LoF variants using the ``--lof`` option.

.. code-block:: bash

	$ gemini pathways --lof -v 68 example.db
	chrom	start	end	ref	alt	impact	sample	genotype	gene	transcript	pathway
	chr10	126678091	126678092	G	A	stop_gain	M128215	G/A	CTBP2	ENST00000531469	hsa05220:Chronic_myeloid_leukemia,hsa04310:Wnt_signaling_pathway,hsa04330:Notch_signaling_pathway,hsa05200:Pathways_in_cancer



===========================================================================
``interactions``: Find genes among variants that are interacting partners.
===========================================================================
Integrating the knowledge of the known protein-protein interactions would be
useful in explaining variation data. Meaning to say that a damaging variant
in an interacting partner of a  potential protein may be equally interesting
as the protein itself. We have used the HPRD binary interaction data to build
a p-p network graph which can be explored by Gemini.


Examples:

.. code-block:: bash

	$ gemini interactions -g CTBP2 -r 3 example.db
	sample	gene	order_of_interaction	interacting_gene
	M128215	CTBP2	0_order:	CTBP2
	M128215	CTBP2	1_order:	RAI2
	M128215	CTBP2	2_order:	RB1
	M128215	CTBP2	3_order:	TGM2,NOTCH2NL

Return CTBP2 (-g) interacting gene variants till the third order (-r)

---------------------
``lof_interactions``
---------------------
Use this option to restrict your analysis to only LoF variants.

.. code-block:: bash

	$ gemini lof_interactions -r 3 example.db
	sample	lof_gene	order_of_interaction	interacting_gene
	M128215	TGM2	1_order:	RB1
	M128215	TGM2	2_order:	none
	M128215	TGM2	3_order:	NOTCH2NL,CTBP2


Meaning to say return all LoF gene TGM2 (in sample M128215) interacting
partners to a 3rd order of interaction.


---------------------
``--var``
---------------------

An extended variant information (chrom, start, end etc.) for the interacting gene
may be achieved with the --var option for both the ``interactions`` and the
``lof_interactions``

.. code-block:: bash

	$ gemini interactions -g CTBP2 -r 3 --var example.db
	sample	gene	order_of_interaction	interacting_gene	var_id	chrom	start	end	impact	biotype	in_dbsnp	clinvar_sig	clinvar_disease_name	aaf_1kg_all	aaf_esp_all
	M128215	CTBP2	0	CTBP2	5	chr10	126678091	126678092	stop_gain	protein_coding	1	None	None	None	None
	M128215	CTBP2	1	RAI2	9	chrX	17819376	17819377	non_syn_coding	protein_coding	1	None	None	1	0.000473
	M128215	CTBP2	2	RB1	7	chr13	48873834	48873835	upstream	protein_coding	1	None	None	0.94	None
	M128215	CTBP2	3	NOTCH2NL	1	chr1	145273344	145273345	non_syn_coding	protein_coding	1	None	None	None	None
	M128215	CTBP2	3	TGM2	8	chr20	36779423	36779424	stop_gain	protein_coding	0	None	None	None	None

.. code-block:: bash

	$ gemini lof_interactions -r 3 --var example.db
	sample	lof_gene	order_of_interaction	interacting_gene	var_id	chrom	start	end	impact	biotype	in_dbsnp	clinvar_sig	clinvar_disease_name	aaf_1kg_all	aaf_esp_all
	M128215	TGM2	1	RB1	7	chr13	48873834	48873835	upstream	protein_coding	1	None	None	0.94	None
	M128215	TGM2	3	NOTCH2NL	1	chr1	145273344	145273345	non_syn_coding	protein_coding	1	None	None	None	None
	M128215	TGM2	3	CTBP2	5	chr10	126678091	126678092	stop_gain	protein_coding	1	None	None	None	None


===================================================================================
``lof_sieve``: Filter LoF variants by transcript position and type
===================================================================================
Not all candidate LoF variants are created equal. For e.g, a nonsense (stop gain)
variant impacting the first 5% of a polypeptide is far more likely to be deleterious
than one affecting the last 5%. Assuming you've annotated your VCF with snpEff v3.0+,
the lof_sieve tool reports the fractional position (e.g. 0.05 for the first 5%) of
the mutation in the amino acid sequence. In addition, it also reports the predicted
function of the transcript so that one can segregate candidate LoF variants that
affect protein_coding transcripts from processed RNA, etc.


.. code-block:: bash

	$ gemini lof_sieve chr22.low.exome.snpeff.100samples.vcf.db
	chrom   start   end ref alt highest_impact  aa_change   var_trans_pos   trans_aa_length var_trans_pct   sample  genotype    gene    transcript  trans_type
	chr22   17072346    17072347    C   T   stop_gain   W365*   365 557 0.655296229803  NA19327 C|T CCT8L2  ENST00000359963 protein_coding
	chr22   17072346    17072347    C   T   stop_gain   W365*   365 557 0.655296229803  NA19375 T|C CCT8L2  ENST00000359963 protein_coding
	chr22   17129539    17129540    C   T   splice_donor    None    None    None    None    NA18964 T|C TPTEP1  ENST00000383140 lincRNA
	chr22   17129539    17129540    C   T   splice_donor    None    None    None    None    NA19675 T|C TPTEP1  ENST00000383140 lincRNA


===========================================================
``annotate``: adding your own custom annotations
===========================================================
It is inevitable that researchers will want to enhance the gemini framework with
their own, custom annotations. ``gemini`` provides a sub-command called
``annotate`` for exactly this purpose. As long as you provide a ``tabix``'ed
annotation file in BED format, the ``annotate`` tool will, for each
variant in the variants table, screen for overlaps in your annotation file and
update a one or more new column in the variants table that you may specify on the command
line. This is best illustrated by example.

Let's assume you have already created a gemini database of a VCF file using
the ``load`` module.

.. code-block:: bash

    $ gemini load -v my.vcf -t snpEff my.db

Now, let's imagine you have an annotated file in BED format (``important.bed``)
that describes regions of the genome that are particularly relevant to your
lab's research. You would like to annotate in the gemini database which variants
overlap these crucial regions. We want to store this knowledge in a new column
in the ``variants`` table called ``important_variant`` that tracks whether a given
variant overlapped (1) or did not overlap (0) intervals in your annotation file.

To do this, you must first TABIX your BED file:

.. code-block:: bash

    $ bgzip important.bed
    $ tabix -p bed important.bed.gz


------------------------------------------------------
``-a boolean`` Did a variant overlap a region or not?
------------------------------------------------------

.. note::

    Formerly, the ``-a`` option was the ``-t`` option.


Now, you can use this TABIX'ed file to annotate which variants overlap your
important regions.  In the example below, the results will be stored in a new
column called "important".  The ``-t boolean`` option says that you just want to
track whether (1) or not (0) the variant overlapped one or more of your regions.

.. code-block:: bash

    $ gemini annotate -f important.bed.gz -c important -a boolean my.db

Since a new columns has been created in the database, we can now directly query
the new column.  In the example results below, the first and third variants
overlapped a crucial region while the second did not.

.. code-block:: bash

    $ gemini query \
        -q "select chrom, start, end, variant_id, important from variants" \
        my.db \
        | head -3
    chr22   100    101    1   1
    chr22   200    201    2   0
    chr22   300    500    3   1


-----------------------------------------------------
``-a count`` How many regions did a variant overlap?
-----------------------------------------------------
Instead of a simple yes or no, we can use the ``-t count`` option to *count*
how many important regions a variant overlapped.  It turns out that the 3rd
variant actually overlapped two important regions.

.. code-block:: bash

    $ gemini annotate -f important.bed.gz -c important -a count my.db

    $ gemini query \
        -q "select chrom, start, end, variant_id, crucial from variants" \
        my.db \
        | head -3
    chr22   100    101    1   1
    chr22   200    201    2   0
    chr22   300    500    3   2


-------------------------------------------------------
``-a extract`` Extract specific values from a BED file
-------------------------------------------------------
Lastly, we may also extract values from specific fields in a BED
file and populate one or more new columns in the database based on
overlaps with the annotation file and the values of the fields therein.
To do this, we use the ``-a extract`` option.

This is best described with an example.  To set this up, let's imagine
that we have a VCF file from a different experiment and we want to annotate
the variants in our GEMINI database with the allele frequency and depth
tags from the INFO fields for the same variants in this other VCF file.

First, since the ``annotate`` tool only supports BED files, we must use the
excellent `vcftools <http://vcftools.sourceforge.net/>`_ package to extract the allele frequency (AF) and
depth (DP) tags from the VCF file.

.. code-block:: bash

    # this will create a new file called other.INFO
    $ vcftools --vcf other.vcf --get-INFO AF --get-INFO DP --out other

    # peek at the output
    $ head -6 other.INFO
    CHROM   POS REF ALT AF  DP
    chr10   1142208 T   C   1.00    122
    chr10   48003992    C   T   0.50    165
    chr10   48004992    C   T   0.50    165
    chr10   135336656   G   A   1.00    2
    chr10   135369532   T   C   0.25    239

    # create a BED file from the output of VCFTOOLs.
    $ awk -v OFS="\t" '{if (NR>1) {print $1,$2-1,$2,$5,$6}}' other.INFO > other.bed

    # peek at the output
    $ head -5 other.bed
    chr10   1142207 1142208 1.00    122
    chr10   48003991    48003992    0.50    165
    chr10   48004991    48004992    0.50    165
    chr10   135336655   135336656   1.00    2
    chr10   135369531   135369532   0.25    239

    # bgzip and tabix for use with the annotate tool.
    $ bgzip other.bed
    $ tabix -p bed other.bed.gz

Now that we have a proper TABIX'ed BED file, we can use the ``-a extract`` option to populate new
columns in the GEMINI database.  In order to do so, we must specify:

    1. the name of the column we want to add (``-c``)

    2. its type (e.g., text, int, float,)  (``-t``)

    3. the column in the BED file that we should use to extract data with which to populate the new column (``-e``)

    4. what operation should be used to summarize the data in the event of multiple overlaps in the annotation file  (``-o``)

For example, let's imagine we want to create a new column called "other_allele_freq" using the
AF column (that is, the 4th column) in our BED file to populate it.

.. code-block:: bash

    $ gemini annotate -f other.bed.gz \
                      -a extract \
                      -c other_allele_freq \
                      -t float \
                      -e 4 \
                      -o mean \
                      my.db

This create a new column in ``my.db`` called ``other_allele_freq`` and this
new column will be a FLOAT.  In the event of multiple records in the BED
file overlapping a variant in the database, the average (mean) of the allele
frequencies values from the BED file will be used.

At this point, one can query the database based on the values of the
new ``other_allele_freq`` column:

.. code-block:: bash

    $ gemini query -q "select * from variants where other_allele_freq < 0.01" my.db


-------------------------------------------------------------------
``-t TYPE`` Specifying the column type(s) when using ``-a extract``
-------------------------------------------------------------------

The ``annotate`` tool will create three different types of columns via the ``-t`` option:

    1. Floating point columns for annotations with decimal precision as above (``-t float``)
    2. Integer columns for integral annotations (``-t integer``)
    3. Text columns for string columns such as "valid", "yes", etc. (``-t text``)

.. note::

    The ``-t`` option is only valid when using the ``-a extract`` option.

----------------------------------------------------------------------------
``-o OPERATION`` Specifying the summary operations when using ``-a extract``
----------------------------------------------------------------------------

In the event of multiple overlaps between a variant and records in the annotation
file, the ``annotate`` tool can summarize the values observed with multiple options:

    1. ``-o mean``.  Compute the average of the values.  **They must be numeric**.
    2. ``-o median``. Compute the median of the values.  **They must be numeric**.
    3. ``-o mix``. Compute the minimum of the values.  **They must be numeric**.
    4. ``-o max``. Compute the maximum of the values.  **They must be numeric**.
    5. ``-o mode``. Compute the maximum of the values.  **They must be numeric**.
    6. ``-o first``. Use the value from the **first** record in the annotation file.
    7. ``-o last``. Use the value from the **last** record in the annotation file.
    8. ``-o list``. Create a comma-separated list of the observed values.  **-t must be text**
    9. ``-o uniq_list``. Create a comma-separated list of the **distinct** (i.e., non-redundant) observed values.  **-t must be text**

.. note::

    The ``-o`` option is only valid when using the ``-a extract`` option.


-------------------------------------------------------------------
Extracting and populating multiple columns at once.
-------------------------------------------------------------------
One can also extract and populate multiple columns at once by providing
comma-separated lists (no spaces) of column names (``-c``), types (``-t``), numbers (``-e``),
and summary operations (``-o``).  For example, recall that in the VCF example above,
we created a TABIX'ed BED file containg the allele frequency and depth values from
the INFO field as the 4th and 5th columns in the BED, respectively.

Instead of running the ``annotate`` tool twice (once for eaxh column), we can
run the tool once and load both columns in the same run.  For example:

.. code-block:: bash

    $ gemini annotate -f other.bed.gz \
                      -a extract \
                      -c other_allele_freq,other_depth \
                      -t float,integer \
                      -e 4,5 \
                      -o mean,max \
                      my.db

We can then use each of the new columns to filter variants with a GEMINI query:

.. code-block:: bash

    $ gemini query -q "select * from variants \
                       where other_allele_freq < 0.01 \
                       and other_depth > 100" my.db


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

---------------------
``--columns``
---------------------

By default, this tool reports all columns in the ``variants`` table. One may
choose to report only a subset of the columns using the ``--columns`` option.  For
example, to report just the ``gene, chrom, start, end, ref, alt, impact``, and ``impact_severity`` columns, one
would use the following:

.. code-block:: bash

    $ gemini region --gene DHODH \
                    --columns "chrom, start, end, ref, alt, gene, impact" \
                    my.db

    chr16   72057281    72057282    A   G   DHODH   intron
    chr16   72057434    72057435    C   T   DHODH   non_syn_coding
    chr16   72059268    72059269    T   C   DHODH   downstream

---------------------
``--filter``
---------------------

By default, this tool will report all variants regardless of their putative
functional impact.  In order to apply additional constraints on the variants
returned, one can use the ``--filter`` option.  Using SQL syntax, conditions
applied with the ``--filter option become WHERE clauses in the query issued to
the GEMINI database.  For example, if we wanted to restrict candidate variants
to solely those with a HIGH predicted functional consequence, we could use the
following:

.. code-block:: bash

    $ gemini region --gene DHODH \
                    --columns "chrom, start, end, ref, alt, gene, impact" \
                    --filter "alt='G'"
                    my.db

    chr16   72057281    72057282    A   G   DHODH   intron

---------------------
``--json``
---------------------
Reporting query output in JSON format may enable
HTML/Javascript apps to query GEMINI and retrieve
the output in a format that is amenable to web development protocols.

To report in JSON format, use the ``--json`` option. For example:

.. code-block:: bash

    $ gemini region --gene DHODH \
                    --columns "chrom, start, end, ref, alt, gene, impact" \
                    --filter "alt='G'"
                    --json
                    my.db

    {"chrom": "chr16", "start": 72057281, "end": 72057282, "ref": "A", "alt": "G", "gene": "DHODH"}



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


Compute the max value for HWE statistic for all variants in a window of size
10kb

.. code-block:: bash

	$ gemini windower  -w 10000 -t hwe -o max my.db


===========================================================================
``stats``: Compute useful variant statistics.
===========================================================================
The ``stats`` tool computes some useful variant statistics like


Compute the transition and transversion ratios for the snps

.. code-block:: bash

	$ gemini stats --tstv my.db
	ts	tv	ts/tv
	4	5	0.8



---------------------
``--tstv-coding``
---------------------
Compute the transition/transversion ratios for the snps in the coding
regions.

----------------------
``--tstv-noncoding``
----------------------
Compute the transition/transversion ratios for the snps in the non-coding
regions.


Compute the type and count of the snps.

.. code-block:: bash

	$ gemini stats --snp-counts my.db
	type	count
	A->G	2
	C->T	1
	G->A	1


Calculate the site frequency spectrum of the variants.

.. code-block:: bash

	$ gemini stats --sfs my.db
	aaf	count
	0.125	2
	0.375	1


Compute the pair-wise genetic distance between each sample

.. code-block:: bash

	$ gemini stats --mds my.db
	sample1	sample2	distance
	M10500	M10500	0.0
	M10475	M10478	1.25
	M10500	M10475	2.0
	M10500	M10478	0.5714



Return a count of the types of genotypes per sample

.. code-block:: bash

	$ gemini stats --gts-by-sample my.db
	sample	num_hom_ref	num_het	num_hom_alt	num_unknown	total
	M10475	4	1	3	1	9
	M10478	2	2	4	1	9



Return the total variants per sample (sum of homozygous
and heterozygous variants)

.. code-block:: bash

	$ gemini stats --vars-by-sample my.db
	sample	total
	M10475	4
	M10478	6


----------------------
``--summarize``
----------------------

If none of these tools are exactly what you want, you can summarize the variants
per sample of an arbitrary query using the --summarize flag. For example, if you
wanted to know, for each sample, how many variants are on chromosome 1 that are also
in dbSNP:

.. code-block:: bash

   	$ gemini stats --summarize "select * from variants where in_dbsnp=1 and chrom='chr1'" my.db
	sample	total	num_het	num_hom_alt
	M10475	1	1	0
	M128215	1	1	0
	M10478	2	2	0
	M10500	2	1	1

===============================================================
``burden``: perform sample-wise gene-level burden calculations
===============================================================
The ``burden`` tool provides a set of utilities to perform burden
summaries on a per-gene, per sample basis. By default, it outputs
a table of gene-wise counts of all high impact variants in coding regions for
each sample:

.. code-block:: bash

	$ gemini burden test.burden.db
	gene	M10475	M10478	M10500	M128215
	WDR37	2	2	2	2
	CTBP2	0	0	0	1
	DHODH	1	0	0	0

----------------------
``--nonsynonymous``
----------------------
If you want to be a little bit less restrictive, you can include all
non-synonymous variants instead:

.. code-block:: bash

   	$ gemini burden --nonsynonymous test.burden.db
	gene	M10475	M10478	M10500	M128215
	SYCE1	0	1	1	0
	WDR37	2	2	2	2
	CTBP2	0	0	0	1
	ASAH2C	2	1	1	0
	DHODH	1	0	0	0

----------------------
``--calpha``
----------------------
If your database has been loaded with a PED file describing case and
control samples, you can calculate the
`c-alpha <http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1001322>`_
statistic for cases vs. control:

.. code-block:: bash

   	$ gemini burden --calpha test.burden.db
	gene	T	c	Z	p_value
	SYCE1	-0.5	0.25	-1.0	0.841344746069
	WDR37	-1.0	1.5	-0.816496580928	0.792891910879
	CTBP2	0.0	0.0	nan	nan
	ASAH2C	-0.5	0.75	-0.57735026919	0.718148569175
	DHODH	0.0	0.0	nan	nan

To calculate the P-value using a permutation test, use the ``--permutations`` option,
specifying the number of permutations of the case/control labels you want to use.

------------------------------------------------
``--min-aaf`` and ``--max-aaf`` for ``--calpha``
------------------------------------------------
By default, all variants affecting a given gene will be included in the
C-alpha computation.  However, one may establish alternate allele frequency
boundaries for the variants included using the ``--min-aaf`` and
``--max-aaf`` options.

.. code-block:: bash

   	$ gemini burden --calpha test.burden.db --min-aaf 0.0 --max-aaf 0.01


---------------------------------------------
``--cases`` and ``--controls for ``--calpha``
---------------------------------------------

If you do not have a PED file loaded, or your PED file does not follow the
standard `PED phenotype encoding format <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml>`_
you can still perform the c-alpha test, but you have to specify which samples
are the control samples and which are the case samples:

.. code-block:: bash

	$ gemini burden --controls M10475 M10478 --cases M10500 M128215 --calpha test.burden.db
	gene	T	c	Z	p_value
	SYCE1	-0.5	0.25	-1.0	0.841344746069
	WDR37	-1.0	1.5	-0.816496580928	0.792891910879
	CTBP2	0.0	0.0	nan	nan
	ASAH2C	-0.5	0.75	-0.57735026919	0.718148569175
	DHODH	0.0	0.0	nan	nan

---------------------------------------------
``--nonsynonymous`` ``--calpha``
---------------------------------------------
If you would rather consider all nonsynonymous variants for the C-alpha test rather
than just the medium and high impact variants, add the ``--nonsynonymous`` flag.


===========================================================================
``ROH``: Identifying runs of homozygosity
===========================================================================
Runs of homozygosity are long stretches of homozygous genotypes that reflect
segments shared identically by descent and are a result of consanguinity or
natural selection. Consanguinity elevates the occurrence of rare recessive 
diseases (e.g. cystic fibrosis) that represent homozygotes for strongly deleterious 
mutations. Hence, the identification of these runs holds medical value. 

The 'roh' tool in GEMINI returns runs of homozygosity identified in whole genome data. 
The tool basically looks at every homozygous position on the chromosome as a possible
start site for the run and looks for those that could give rise to a potentially long 
stretch of homozygous genotypes. 

For e.g. for the given example allowing ``1 HET`` genotype (h) and ``2 UKW`` genotypes (u) 
the possible roh runs (H) would be:


.. code-block:: bash

	genotype_run = H H H H h H H H H u H H H H H u H H H H H H H h H H H H H h H H H H H
	roh_run1     = H H H H h H H H H u H H H H H u H H H H H H H
	roh_run2     =           H H H H u H H H H H u H H H H H H H h H H H H H
	roh_run3     =                     H H H H H u H H H H H H H h H H H H H
	roh_run4     =                                 H H H H H H H h H H H H H

roh returned for --min-snps = 20 would be:

.. code-block:: bash
	
	roh_run1     = H H H H h H H H H u H H H H H u H H H H H H H
	roh_run2     =           H H H H u H H H H H u H H H H H H H h H H H H H

As you can see, the immediate homozygous position right of a break (h or u) would be the possible 
start of a new roh run and genotypes to the left of a break are pruned since they cannot 
be part of a longer run than we have seen before.



Return ``roh`` with minimum of 50 snps, a minimum run length of 1 mb and a minimum sample depth of 20 
for sample S138 (with default values for allowed number of HETS, UNKS and total depth). 

.. code-block:: bash

	$ gemini roh --min-snps 50 \
	           --min-gt-depth 20 \
			   --min-size 1000000 \
			   -s S138 \
			   roh_run.db
	chrom	start	end	sample	num_of_snps	density_per_kb	run_length_in_bp
	chr2 233336080 234631638 S138 2583 1.9953 1295558
	chr2	238341281	239522281	S138	2899	2.4555	1181000


===========================================================================
``set_somatic``: Flag somatic variants
===========================================================================
Somatic mutations in a tumor-normal pair are variants that are present in
the tumor but not in the normal sample.

.. note::

    1. This tool requires that you specify the sample layout via a PED file
    when loading your VCF into GEMINI via:

    ``gemini load -v my.vcf -p my.ped my.db``


`Example PED file format for GEMINI`

.. code-block:: bash

	#Family_ID	Individual_ID	Paternal_ID	Maternal_ID	Sex	Phenotype	Ethnicity
	1       Normal  -9      -9      0       1       -9
	1       Tumor   -9      -9      0       2       -9


---------------------
``default behavior``
---------------------
By default, ``set_somatic`` simply marks variants that are genotyped as
homozygous reference in the normal sample and non-reference in the tumor. 
More stringent somatic filtering criteria are available through tunable
command line parameters.

.. code-block:: bash

	$ gemini set_somatic \
            --min-depth 30 \
            --min-qual 20 \
            --min-somatic-score 18 \
            --min-tumor-depth 10 \
            --min-norm-depth 10 \
            tumor_normal.db
        tum_name	tum_gt	tum_alt_freq	tum_alt_depth	tum_depth	nrm_name	nrm_gt	nrm_alt_freq	nrm_alt_depth	nrm_depth	chrom	start	end	ref	alt	gene
        tumor	GAAAAAAAAAAAAAGGTGAAAATT/GAAAAAAAAAAAAGGTGAAAATT	0.217391304348	5	23	normal	GAAAAAAAAAAAAAGGTGAAAATT/GAAAAAAAAAAAAAGGTGAAAATT	0.0	0	25	chrX	132838304	132838328	GAAAAAAAAAAAAAGGTGAAAATT	GAAAAAAAAAAAAGGTGAAAATT	GPC3
        tumor	CTGCTATTTTG/CG	0.22	11	50	normal	CTGCTATTTTG/CTGCTATTTTG	0.0	0	70	chr17	59861630	59861641	CTGCTATTTTG	CG	BRIP1
        tumor	C/A	0.555555555556	10	18	normal	C/C	0.0	0	17	chr17	7578460	7578461	C	A	TP53
        tumor	C/T	0.1875	12	64	normal	C/C	0.0	0	30	chr2	128046288	128046289	C	T	ERCC3
        Identified and set 4 somatic mutations


---------------------
``--min-depth [None]``
---------------------
The minimum required combined depth for tumor and normal samples.

---------------------
``--min-qual [None]``
---------------------
The minimum required variant quality score.

---------------------
``--min-somatic-score [None]``
---------------------
The minimum required somatic score (SSC). This score is produced by various
somatic variant detection algorithms including SpeedSeq, SomaticSniper,
and VarScan 2.

---------------------
``--max-norm-alt-freq [None]``
---------------------
The maximum frequency of the alternate allele allowed in the normal sample.

---------------------
``--max-norm-alt-count [None]``
---------------------
The maximum count of the alternate allele allowed in the normal sample.

---------------------
``--min-norm-depth [None]``
---------------------
The minimum depth required in the normal sample.

---------------------
``--min-tumor-alt-freq [None]``
---------------------
The minimum frequency of the alternate allele required in the tumor sample.

---------------------
``--min-tumor-alt-count [None]``
---------------------
The minimum count of the alternate allele required in the tumor sample.

---------------------
``--min-tumor-depth [None]``
---------------------
The minimum depth required in the tumor sample.

---------------------
``--chrom [None]``
---------------------
A specific chromosome on which to flag somatic mutations.

---------------------
``--dry-run``
---------------------
Don't set the is_somatic flag, just report what _would_ be set. For testing
purposes.


===========================================================================
``actionable_mutations``: Report actionable somatic mutations and drug-gene interactions
===========================================================================
Actionable mutations are somatic variants in COSMIC cancer census genes with
medium or high impact severity predictions. This tool reports actionable
mutations as well as their known drug interactions (if any) from DGIdb.
Current functionality is only for SNVs and indels.

.. note::

    1. This tool requires somatic variants to have been flagged using
           ``set_somatic``


.. code-block:: bash

	$ gemini actionable_mutations tumor_normal.db
	tum_name	chrom	start	end	ref	alt	gene	impact	is_somatic	in_cosmic_census	dgidb_info
	tumor	chr2	128046288	128046289	C	T	ERCC3	non_syn_coding	1	1	None
	tumor	chr17	7578460	7578461	C	A	TP53	non_syn_coding	1	1	{'searchTerm': 'TP53', 'geneCategories': ['CLINICALLY ACTIONABLE', 'DRUGGABLE GENOME', 'TUMOR SUPPRESSOR', 'TRANSCRIPTION FACTOR COMPLEX', 'DRUG RESISTANCE', 'HISTONE MODIFICATION', 'DNA REPAIR', 'TRANSCRIPTION FACTOR BINDING'], 'geneName': 'TP53', 'geneLongName': 'tumor protein p53', 'interactions': [{'source': 'DrugBank', 'interactionId': '711cbe42-4930-4b46-963e-79ab35bbbd0f', 'interactionType': 'n/a', 'drugName': '1-(9-ETHYL-9H-CARBAZOL-3-YL)-N-METHYLMETHANAMINE'}, {'source': 'PharmGKB', 'interactionId': '8234d9b9-085d-49b1-aac2-cf5375d91477', 'interactionType': 'n/a', 'drugName': 'FLUOROURACIL'}, {'source': 'PharmGKB', 'interactionId': '605d7bca-7ed9-428e-aa7c-f76aafd66b54', 'interactionType': 'n/a', 'drugName': 'PACLITAXEL'}, {'source': 'TTD', 'interactionId': '1fe9db63-3581-435b-b22a-12d45c8c9864', 'interactionType': 'activator', 'drugName': 'CURAXIN CBLC102'}, {'source': 'TALC', 'interactionId': '8f8f6822-cb9e-40aa-8360-5532e059f1e7', 'interactionType': 'vaccine', 'drugName': 'EP-2101'}, {'source': 'TALC', 'interactionId': 'd59e14bc-b9a5-4c9f-a5aa-7ba322f0fa0e', 'interactionType': 'vaccine', 'drugName': 'MUTANT P53 PEPTIDE PULSED DENDRITIC CELL'}, {'source': 'TALC', 'interactionId': '79256b6e-9a16-4fbe-a237-28dbca28bc2a', 'interactionType': 'vaccine', 'drugName': 'AD.P53-DC'}]}
	tumor	chr17	59861630	59861641	CTGCTATTTTG	CG	BRIP1	inframe_codon_loss	1	1	None
	tumor	chrX	132838304	132838328	GAAAAAAAAAAAAAGGTGAAAATT	GAAAAAAAAAAAAGGTGAAAATT	GPC3	splice_region	1	1	None


===========================================================================
``fusions``: Report putative gene fusions
===========================================================================
Report putative somatic gene fusions from structural variants in a tumor-normal
pair. Putative fusions join two genes and preserve transcript strand
orientation.

.. note::

    1. This tool requires somatic variants to have been flagged using
           ``set_somatic``


---------------------
``default behavior``
---------------------
By default, ``fusions`` reports structural variants that are flagged as
somatic, join two different genes, and preserve transcript strand orientation.
These may be further filtered using tunable command line parameters.


.. code-block:: bash

	$ gemini fusions \
	    --min_qual 5 \
	    --in_cosmic_census \
	    tumor_normal.db
	chromA   breakpointA_start  breakpointA_end	chromB	breakpointB_start   breakpointB_end var_id  qual    strandA strandB sv_type geneA   geneB   tool    evidence_type   is_precise  sample
    chr3	176909953	176909982	chr3	178906001	178906030	1233	9.58	-	+	complex	TBL1XR1	PIK3CA	LUMPY	PE	0	tumor


---------------------
``--min_qual [None]``
---------------------
The minimum required variant quality score.

---------------------
``--evidence_type STRING``
---------------------
The required supporting evidence types for the variant from
LUMPY ("PE", "SR", or "PE,SR").

---------------------
``--in_cosmic_census``
---------------------
Require at least one of the affected genes to be in the
COSMIC cancer gene census.


===========================================================================
``db_info``: List the gemini database tables and columns
===========================================================================

Because of the sheer number of annotations that are stored in gemini, there are
admittedly too many columns to remember by rote.  If you can't recall the name of
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
