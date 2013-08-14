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
	NA19675  PKDREJ  chr22,46653547,46653548,C,T,C|T,non_syn_coding,exon_22_46651560_46659219,0.005,1  chr22,46657894,46657895,G,A,A|G,non_syn_coding,exon_22_46651560_46659219,0.005,1


This indicates that sample NA19675 has a candidate compound heterozygote in
PKDREJ.  The two heterozygotes are reported using the following structure:

``chrom,start,end,ref,alt,genotype,impact,exon,AAF,in_dbsnp``

---------------
``--only_lof``
---------------
By default, all coding variants are explored.  However, one may want to
restrict the analysis to LoF variants using the ``--only_lof`` option.

.. code-block:: bash

	$ gemini comp_hets --only_lof chr22.low.exome.snpeff.100samples.vcf.db
	NA19002	GTSE1	chr22,46722400,46722401,G,A,G|A,stop_gain,exon_22,0.005,1	chr22,46704499,46704500,C,A,A|C,stop_gain,exon_22,0.005,0


--------------------
``--only-affected``
--------------------
By default, candidate compound heterozygous variants are reported for all
individuals in the database.  One can restrict the analysis to variants in
only individuals with an affected phenotype using the ``--only-affected`` option.

.. code-block:: bash

	$ gemini comp_hets --only-affected chr22.low.exome.snpeff.100samples.vcf.db


----------------------
``--allow-other-hets``
----------------------
By default, the ``comp_hets`` tool will identify candidate pairs of
heterozygotes where each variant is found in *only one* of the samples in your database.
Depending on the genetic model, this may be too restrictive.  If you'd like to
include variants for which other individuals may also be heterozygous, just use
the ``--allow-other-hets`` option

.. code-block:: bash

	$ gemini comp_hets --allow-other-hets chr22.low.exome.snpeff.100samples.vcf.db
	NA19375  PKDREJ  chr22,46658977,46658978,T,C,T|C,non_syn_coding,exon_22_46651560_46659219,0.25,1  chr22,46655778,46655779,G,C,C|G,non_syn_coding,exon_22_46651560_46659219,0.08,1
	HG01619  PKDREJ  chr22,46658977,46658978,T,C,C|T,non_syn_coding,exon_22_46651560_46659219,0.25,1  chr22,46657307,46657308,T,C,T|C,non_syn_coding,exon_22_46651560_46659219,0.005,1

Here, samples NA19375 and HG01619 are both hets for the same variant (chr22,46658977,46658978)


---------------------
``--ignore-phasing``
---------------------
If your genotypes aren't phased, we can't be certain that two heterozygotes
are on opposite alleles.  However, we can still identify pairs of heterozygotes
that are *candidates* for compound heterozygotes. Just use the
``--ignore-phasing`` option.

.. code-block:: bash

	$ gemini comp_hets --ignore_phasing example.db
	M1047  DHODH  chr16,72048539,72048540,C,T,C/T,non_syn_coding,3/4,0.125,1  chr16,72057434,72057435,C,T,C/T,non_syn_coding,8/9,0.125,1
	M1282  DHODH  chr16,72055099,72055100,C,T,C/T,non_syn_coding,5/9,0.125,0  chr16,72055114,72055116,CT,C,CT/C,frame_shift,5/9,0.125,0



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
family in the database, a all columns in the variants table for mutations that
are not found in the parents yet are observed as heterozygotes in the offspring.
For example:

.. code-block:: bash

    $ gemini de_novo my.db

    family_id	family_members	genotypes	depths	chrom	start	end	variant_id	anno_id	ref	alt	qual	filter	type	sub_type	call_rate	in_dbsnp	rs_ids	in_omim	clinvar_sig	clinvar_disease_name	clinvar_dbsource	clinvar_dbsource_id	clinvar_origin	clinvar_dsdb	clinvar_dsdbid	clinvar_disease_acc	clinvar_in_locus_spec_db	clinvar_on_diag_assay	pfam_domain	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	gerp_bp_score	gerp_element_pval	num_hom_ref	num_het	num_hom_alt	num_unknown	aaf	hwe	inbreeding_coeff	pi	recomb_rate	gene	transcript	is_exonic	is_coding	is_lof	exon	codon_change	aa_change	aa_length	biotype	impact	impact_severity	polyphen_pred	polyphen_score	sift_pred	sift_score	anc_allele	rms_bq	cigar	depth	strand_bias	rms_map_qual	in_hom_run	num_mapq_zero	num_alleles	num_reads_w_dels	haplotype_score	qual_depth	allele_count	allele_bal	in_hm2	in_hm3	is_somatic	in_esp	aaf_esp_ea	aaf_esp_aa	aaf_esp_all	exome_chip	in_1kg	aaf_1kg_amr	aaf_1kg_asn	aaf_1kg_afr	aaf_1kg_eur	aaf_1kg_all	grc	gms_illumina	gms_solid	gms_iontorrent	in_cse	encode_tfbs	encode_dnaseI_cell_count	encode_dnaseI_cell_list	encode_consensus_gm12878	encode_consensus_h1hesc	encode_consensus_helas3	encode_consensus_hepg2	encode_consensus_huvec	encode_consensus_k562	gts	gt_types	gt_phases	gt_depths	gt_ref_depths	gt_alt_depths	gt_quals
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

    family_id	family_members	family_genotypes	depths	chrom	start	end	variant_id	anno_id	ref	alt	qual	filter	type	sub_type	call_rate	in_dbsnp	rs_ids	in_omim	clinvar_sig	clinvar_disease_name	clinvar_dbsource	clinvar_dbsource_id	clinvar_origin	clinvar_dsdb	clinvar_dsdbid	clinvar_disease_acc	clinvar_in_locus_spec_db	clinvar_on_diag_assay	pfam_domain	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	gerp_bp_score	gerp_element_pval	num_hom_ref	num_het	num_hom_alt	num_unknown	aaf	hwe	inbreeding_coeff	pi	recomb_rate	gene	transcript	is_exonic	is_coding	is_lof	exon	codon_change	aa_change	aa_length	biotype	impact	impact_severity	polyphen_pred	polyphen_score	sift_pred	sift_score	anc_allele	rms_bq	cigar	depth	strand_bias	rms_map_qual	in_hom_run	num_mapq_zero	num_alleles	num_reads_w_dels	haplotype_score	qual_depth	allele_count	allele_bal	in_hm2	in_hm3	is_somatic	in_esp	aaf_esp_ea	aaf_esp_aa	aaf_esp_all	exome_chip	in_1kg	aaf_1kg_amr	aaf_1kg_asn	aaf_1kg_afr	aaf_1kg_eur	aaf_1kg_all	grc	gms_illumina	gms_solid	gms_iontorrent	in_cse	encode_tfbs	encode_dnaseI_cell_count	encode_dnaseI_cell_list	encode_consensus_gm12878	encode_consensus_h1hesc	encode_consensus_helas3	encode_consensus_hepg2	encode_consensus_huvec	encode_consensus_k562	gts	gt_types	gt_phases	gt_depths	gt_ref_depths	gt_alt_depths	gt_quals
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

    family_id	family_members	family_genotypes	depths	chrom	start	end	ref	alt
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

    family_id	family_members	family_genotypes	depths	chrom	start	end	ref	alt
    3	243(father; unknown),244(mother; unknown),242(child; affected)	C/C,C/C,C/A	249,243,250	chr1	17729	17730	C	A
    4	254(father; unknown),255(mother; unknown),253(child; affected)	A/A,A/A,A/G	86,146,83	chr1	168097	16809	A	G
    4	254(father; unknown),255(mother; unknown),253(child; affected)	G/G,G/G,G/T	107,182,72	chr1	12854400	12854401	G	T
    3	243(father; unknown),244(mother; unknown),242(child; affected)	A/A,A/A,A/ATGGTGTTG	211,208,208	chr1	12855995	12855996	A	ATGGTGTTG
    ...

============================================================================
``autosomal_recessive``: Find variants meeting an autosomal recessive model.
============================================================================
.. note::

    This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:

    ``gemini load -v my.vcf -p my.ped my.db``

---------------------
``default behavior``
---------------------

Assuming you have defined the familial relationships between samples when loading
your VCF into GEMINI, one can leverage a built-in tool for identifying variants
that meet an autosomal recessive inheritance pattern. The reported variants
will be restricted to those variants having the potential to impact the
function of affecting protein coding transcripts.

.. code-block:: bash

    $ gemini autosomal_recessive my.db

    family_id	family_members	family_genotypes	chrom	start	end	variant_id	anno_id	ref	alt	qual	filter	type	sub_type	call_rate	in_dbsnp	rs_ids	in_omim	clinvar_sig	clinvar_disease_name	clinvar_dbsource	clinvar_dbsource_id	clinvar_origin	clinvar_dsdb	clinvar_dsdbid	clinvar_disease_acc	clinvar_in_locus_spec_db	clinvar_on_diag_assay	pfam_domain	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	gerp_bp_score	gerp_element_pval	num_hom_ref	num_het	num_hom_alt	num_unknown	aaf	hwe	inbreeding_coeff	pi	recomb_rate	gene	transcript	is_exonic	is_coding	is_lof	exon	codon_change	aa_change	aa_length	biotype	impact	impact_severity	polyphen_pred	polyphen_score	sift_pred	sift_score	anc_allele	rms_bq	cigar	depth	strand_bias	rms_map_qual	in_hom_run	num_mapq_zero	num_alleles	num_reads_w_dels	haplotype_score	qual_depth	allele_count	allele_bal	in_hm2	in_hm3	is_somatic	in_esp	aaf_esp_ea	aaf_esp_aa	aaf_esp_all	exome_chip	in_1kg	aaf_1kg_amr	aaf_1kg_asn	aaf_1kg_afr	aaf_1kg_eur	aaf_1kg_all	grc	gms_illumina	gms_solid	gms_iontorrent	in_cse	encode_tfbs	encode_dnaseI_cell_count	encode_dnaseI_cell_list	encode_consensus_gm12878	encode_consensus_h1hesc	encode_consensus_helas3	encode_consensus_hepg2	encode_consensus_huvec	encode_consensus_k562	gts	gt_types	gt_phases	gt_depths	gt_ref_depths	gt_alt_depths	gt_quals
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	CAG/C,CAG/C,C/C	chr1	13655	13658	11	12	CAG	C	1113.97	None	indel	del	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	0	10	2	0	0.583333333333	0.0133475753029	-0.714285714286	0.507246376812	2.981822	DDX11L1	ENST00000518655	0	0	0	3	None	None	None	transcribed_unprocessed_pseudogene	splice_donor	HIGH	None	None	None	None	None	None	None	208	None	17.12	0	1	24	None	19.8498	5.36	14	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	0.0	0.0	21.1	0	None	None	None	R	R	unknown	R	unknown	T	CAG/C,CAG/C,CAG/C,CAG/C,CAG/C,CAG/C,CAG/C,CAG/C,C/C,C/C,CAG/C,CAG/C	1,1,1,1,1,1,1,1,3,3,1,1	False,False,False,False,False,False,False,False,False,False,False,False	9,7,19,13,3,5,4,2,1,1,3,4	14,20,28,29,6,7,6,4,5,6,8,12	13,3,4,10,4,7,4,3,7,2,3,3	99.0,27.88,89.69,99.0,20.9,18.78,77.9,42.91,3.0,3.0,39.85,36.89
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	C/G,C/G,G/G	chr1	726943	726944	283	1	C	G	2492.27	None	snp	tv	1.0	1	rs3131979	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	Satellite_Satellite_(GAATG)n;trf;trf;trf;Satellite_Satellite_(GAATG)n;Satellite_Satellite_(GAATG)n;Satellite_Satellite_(GAATG)n;Satellite_Satellite_(GAATG)n;trf	0	1	0	None	None	0	8	4	0	0.666666666667	0.0832645169833	-0.5	0.463768115942	2.671779	AL669831.1	ENST00000358533	0	0	0	None	None	None	168	protein_coding	downstream	LOW	None	None	None	None	None	None	None	116	None	63.69	1	2	24	0.0	5.4008	21.49	16	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	None	None	None	None	None	None	C/G,C/G,G/G,C/G,G/G,C/G,G/G,C/G,C/G,G/G,C/G,C/G	1,1,3,1,3,1,3,1,1,3,1,1	False,False,False,False,False,False,False,False,False,False,False,False	22,13,9,16,6,9,7,8,4,4,6,12	8,1,1,7,0,7,0,1,3,0,2,7	13,12,8,9,5,1,7,7,1,4,4,5	99.0,0.22,4.39,99.0,15.05,16.57,21.05,11.18,24.09,12.04,19.31,99.0
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	T/C,T/C,C/C	chr1	909308	909309	639	6	T	C	575.77	None	snp	ts	1.0	1	rs3829738	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	0	0	None	7.73127e-14	8	3	1	0	0.208333333333	0.401650457515	0.242105263158	0.344202898551	0.952858	PLEKHN1	ENST00000379407	1	1	0	13	Tcc/Ccc	S476P	576	protein_coding	non_syn_coding	MED	None	None	None	None	None	None	None	274	None	69.97	2	0	24	0.0	0.213	10.28	5	None	None	None	None	1	0.170433	0.24261	0.194885	1	1	0.23	0.27	0.21	0.18	0.22	None	None	None	None	0	None	None	None	R	unknown	T	T	R	R	T/T,T/T,T/T,T/T,T/T,T/T,T/C,T/T,T/T,C/C,T/C,T/C	0,0,0,0,0,0,1,0,0,3,1,1	False,False,False,False,False,False,False,False,False,False,False,False	51,22,36,43,13,19,17,17,17,8,14,17	51,22,36,42,12,19,15,16,17,0,7,9	0,0,0,1,0,0,2,1,0,8,7,8	99.0,57.17,84.27,69.02,33.1,51.15,8.66,2.98,42.12,18.05,99.0,99.0
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	GT/G,GT/G,G/G	chr1	970561	970563	770	2	GT	G	370.77	None	indel	del	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	0	0	None	None	1	10	1	0	0.5	0.0209213346218	-0.666666666667	0.521739130435	0.521736	AGRN	ENST00000379370	0	0	0	2	None	None	2045	protein_coding	intron	LOW	None	None	None	None	None	None	None	83	None	70.33	1	0	24	None	31.3087	4.63	12	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	unknown	CTCF	CTCF	T	unknown	unknown	GT/G,GT/G,GT/G,GT/G,GT/G,GT/G,GT/G,GT/GT,GT/G,G/G,GT/G,GT/G	1,1,1,1,1,1,1,0,1,3,1,1	False,False,False,False,False,False,False,False,False,False,False,False	19,7,8,12,4,6,7,3,6,2,3,6	13,4,5,9,4,5,7,3,4,1,2,5	1,0,3,0,0,1,0,0,2,1,1,0	97.36,28.09,77.28,42.55,16.57,24.71,2.62,8.98,52.45,4.56,19.51,16.39
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
example, to report just the ``chrom, start, end, ref``, and ``alt`` columns, one
would use the following:

.. code-block:: bash

    $ gemini autosomal_recessive --columns "chrom, start, end, ref, alt" my.db

    family_id	family_members	family_genotypes	chrom	start	end	ref	alt
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	CAG/C,CAG/C,C/C	chr1	13655	13658	CAG	C
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	C/G,C/G,G/G	chr1	726943	726944	C	G
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	T/C,T/C,C/C	chr1	909308	909309	T	C
    4	254(father; unaffected),255(mother; unaffected),253(child; affected)	GT/G,GT/G,G/G	chr1	970561	970563	GT	G
    ...

.. note::

    The output will always start with the family ID, the family members, and the
    observed genotypes for the family members.


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
          --columns "chrom, start, end, ref, alt" \
          --filter "impact_severity = 'HIGH'" \
          my.db

    family_id	family_members	family_genotypes	impact_severity	chrom	start	end	ref	alt
    4	SMS254(father; unaffected),SMS255(father; unaffected),SMS253(child; affected)	CAG/C,CAG/C,C/C	chr1	13655	13658	CAG	C
    4	SMS254(father; unaffected),SMS255(father; unaffected),SMS253(child; affected)	A/T,A/T,T/T	chr1	5935161	5935162	A	T
    4	SMS254(father; unaffected),SMS255(father; unaffected),SMS253(child; affected)	C/CGT,C/CGT,CGT/CGT	chr1	20020993	20020994	C	CGT
    4	SMS254(father; unaffected),SMS255(father; unaffected),SMS253(child; affected)	G/GTG,G/GTG,GTG/GTG	chr1	20020994	20020995	G	GTG
    ...





===========================================================================
``autosomal_dominant``: Find variants meeting an autosomal dominant model.
===========================================================================
.. note::

    1. This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:

    ``gemini load -v my.vcf -p my.ped my.db``

    2. If neither parent are known to be affected, this tool will report any
       variant where one and only of the parents is heterozygous and the affected
       child is also heterozygous.  If one and only one of the parents is affected,
       the tool will report variants where both the affected child and the affected
       parent are heterozygous.  If both parents are known to be affected, the
       tool will report nothing for that family.

---------------------
``default behavior``
---------------------

Assuming you have defined the familial relationships between samples when loading
your VCF into GEMINI, one can leverage a built-in tool for identifying variants
that meet an autosomal dominant inheritance pattern. The reported variants
will be restricted to those variants having the potential to impact the
function of affecting protein coding transcripts.

.. code-block:: bash

    $ gemini autosomal_dominant my.db | head

    family_id	family_members	family_genotypes	chrom	start	end	variant_id	anno_id	ref	alt	qual	filter	type	sub_type	call_rate	in_dbsnp	rs_ids	in_omim	clinvar_sig	clinvar_disease_name	clinvar_dbsource	clinvar_dbsource_id	clinvar_origin	clinvar_dsdb	clinvar_dsdbid	clinvar_disease_acc	clinvar_in_locus_spec_db	clinvar_on_diag_assay	pfam_domain	cyto_band	rmsk	in_cpg_island	in_segdup	is_conserved	gerp_bp_score	gerp_element_pval	num_hom_ref	num_het	num_hom_alt	num_unknown	aaf	hwe	inbreeding_coeff	pi	recomb_rate	gene	transcript	is_exonic	is_coding	is_lof	exon	codon_change	aa_change	aa_length	biotype	impact	impact_severity	polyphen_pred	polyphen_score	sift_pred	sift_score	anc_allele	rms_bq	cigar	depth	strand_bias	rms_map_qual	in_hom_run	num_mapq_zero	num_alleles	num_reads_w_dels	haplotype_score	qual_depth	allele_count	allele_bal	in_hm2	in_hm3	is_somatic	in_esp	aaf_esp_ea	aaf_esp_aa	aaf_esp_all	exome_chip	in_1kg	aaf_1kg_amr	aaf_1kg_asn	aaf_1kg_afr	aaf_1kg_eur	aaf_1kg_all	grc	gms_illumina	gms_solid	gms_iontorrent	in_cse	encode_tfbs	encode_dnaseI_cell_count	encode_dnaseI_cell_list	encode_consensus_gm12878	encode_consensus_h1hesc	encode_consensus_helas3	encode_consensus_hepg2	encode_consensus_huvec	encode_consensus_k562	gts	gt_types	gt_phases	gt_depths	gt_ref_depths	gt_alt_depths	gt_quals
    4	SMS254(father; unknown),SMS255(mother; unknown),SMS253(child; affected)	G/G,G/C,G/C	chr1	13272	13273	7	1	G	C	1647.92	None	snp	tv	1.0	0	None	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	4.21522e-07	7	5	0	0	0.208333333333	0.36197632685	-0.263157894737	0.344202898551	2.981822	WASH7	ENST00000423562	0	0	0	None	None	None	None	unprocessed_pseudogene	downstream	LOW	None	None	None	None	None	None	None	649	None	24.18	0	18	24	0.0	1.470	7.26	5	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	44.7	45.2	73.8	0	None	None	None	R	R	unknown	R	unknown	T	G/C,G/G,G/G,G/G,G/C,G/C,G/G,G/G,G/G,G/C,G/G,G/C	1,0,0,0,1,1,0,0,0,1,0,1	False,False,False,False,False,False,False,False,False,False,False,False	74,63,135,83,39,52,31,43,22,23,45,39	46,62,135,81,29,35,31,42,22,15,43,25	28,1,0,2,10,17,0,0,0,8,1,14	99.0,87.24,99.0,75.19,99.0,99.0,42.12,36.1,12.03,99.0,29.72,99.0
    2	SMS230(father; unaffected),SMS231(mother; affected),SMS193(child; affected)	C/C,C/T,C/T	chr1	13301	13302	8	1	C	T	39.25	None	snp	ts	1.0	1	rs180734498	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	4.21522e-07	5	7	0	0	0.291666666667	0.15375441881	-0.411764705882	0.43115942029	2.981822	WASH7P	ENST00000423562	0	0	0	None	None	None	None	unprocessed_pseudogene	downstream	LOW	None	None	None	None	None	None	None	487	None	30.03	0	13	24	0.0	0.8688	0.15	7	None	None	None	None	0	None	None	None	0	1	0.08	0.02	0.21	0.14	0.11	None	None	None	None	0	None	None	None	R	R	unknown	R	unknown	T	C/T,C/T,C/C,C/T,C/T,C/T,C/C,C/C,C/C,C/T,C/C,C/T	1,1,0,1,1,1,0,0,0,1,0,1	False,False,False,False,False,False,False,False,False,False,False,False	56,49,119,55,23,43,23,27,10,18,38,26	36,36,107,28,13,30,16,12,1,12,32,15	20,13,12,27,10,13,7,15,9,6,6,11	4.6,11.62,99.0,0.68,19.06,8.5,33.07,30.07,3.01,22.09,72.16,36.78
    1	SMS238(father; affected),SMS239(father; unaffected),SMS173(child; affected)	G/A,G/G,G/A	chr1	14975	14976	29	1	G	A	5222.86	None	snp	ts	1.0	1	rs71252251	None	None	None	None	None	None	None	None	None	None	None	None	chr1p36.33	None	0	1	0	None	None	4	8	0	0	0.333333333333	0.0832645169833	-0.5	0.463768115942	2.981822	DDX11L1	ENST00000456328	0	0	0	None	None	None	None	processed_transcript	downstream	LOW	None	None	None	None	None	None	None	2993	None	34.6	0	102	24	0.0	3.7036	2.62	8	None	None	None	None	0	None	None	None	0	0	None	None	None	None	None	None	None	None	None	0	None	None	None	R	R	R	R	R	T	G/A,G/A,G/G,G/A,G/A,G/G,G/A,G/G,G/A,G/A,G/G,G/A	1,1,0,1,1,0,1,0,1,1,0,1	False,False,False,False,False,False,False,False,False,False,False,False	250,250,250,250,250,250,250,250,248,245,250,250	211,220,245,204,162,250,196,248,209,216,249,210	38,29,5,46,88,0,53,2,39,28,1,40	99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0,99.0


---------------------
``--columns``
---------------------

By default, this tool reports all columns in the ``variants`` table. One may
choose to report only a subset of the columns using the ``--columns`` option.  For
example, to report just the ``chrom, start, end, ref``, and ``alt`` columns, one
would use the following:

.. code-block:: bash

    $ gemini autosomal_dominant --columns "chrom, start, end, ref, alt" my.db

    family_id	family_members	family_genotypes	chrom	start	end	ref	alt
    4	SMS254(father; unknown),SMS255(mother; unknown),SMS253(child; affected)	G/G,G/C,G/C	chr1	13272	13273	G	C
    2	SMS230(father; unaffected),SMS231(mother; affected),SMS193(child; affected)	C/C,C/T,C/T	chr1	13301	13302	C	T
    1	SMS238(father; affected),SMS239(father; unaffected),SMS173(child; affected)	G/A,G/G,G/A	chr1	14975	14976	G	A
    ...

.. note::

    The output will always start with the family ID, the family members, and the
    observed genotypes for the family members.


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
          --columns "chrom, start, end, ref, alt" \
          --filter "impact_severity = 'HIGH'" \
          my.db

    family_id	family_members	family_genotypes	impact_severity	chrom	start	end	ref	alt
    4	SMS254(father; unknown),SMS255(father; unknown),SMS253(child; affected)	TTCT/T,TTCT/TTCT,TTCT/T	chr1	17362	17366	TTCT	T
    1	SMS238(father; affected),SMS239(father; unaffected),SMS173(child; affected)	G/A,G/G,G/A	chr1	13580	135804	G	A
    1	SMS238(father; affected),SMS239(father; unaffected),SMS173(child; affected)	G/C,C/C,G/C	chr1	99858	998582	G	C
    2	SMS230(father; unaffected),SMS231(father; affected),SMS193(child; affected)	C/C,G/C,G/C	chr1	99858    ...
    ...


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

.. note::

	The ``--calpha`` option merely computes an asymptotic P-value that
	assumes a normal distribution.  It does not yet perform permutation
	tests to yield appropriate p-value distributions even in the presence of LD between variants.

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
	DHODH	0.0	0.0	nan	nan" > exp


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
