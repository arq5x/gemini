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

    This tool requires that you identify familial relationships via a PED file
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
family in the database, a list of mutations that are not found in the parents yet
are observed as heterozygotes in the offspring. For example:

.. code-block:: bash

	$ gemini de_novo my.db

	family_id	chrom	start	end	ref	alt	gene	impact	impact_severity	in_dbsnp	rs_ids		aaf_1kg_all	aaf_esp_all	clinvar_sig	clinvar_disease_name	clinvar_dbsource	sample1(father)	sample2(mother)	sample3(child; affected)	sample1(depth)	sample2(depth)	sample3(depth)
	1	chr1	17197609	17197610	G	A	BX284668.1	non_syn_coding	MED	1	rs200754171		None	None	None	None	None	G/G	G/G	G/A	104	168	244
	1	chr1	196763706	196763707	T	C	CFHR3	splice_acceptor	HIGH	1	rs481759		None	None	None	None	None	T/T	T/T	T/C	26	28	34
	1	chr1	248813541	248813542	G	A	OR2T27	non_syn_coding	MED	1	rs77685347	0.	17	0.180025	None	None	None	G/G	G/G	G/A	21	38	68
	1	chr2	90060872	90060873	A	T	AC009958.1	non_syn_coding	MED	1	rs202041075		None	None	None	None	None	A/A	A/A	A/T	90	238	234
	1	chr3	195505789	195505790	G	C	MUC4	non_syn_coding	MED	1	rs11928301		None	None	None	None	None	G/G	G/G	G/C	250	247	248
	...

---------------------
``-d``
---------------------

Unfortunately, inherited variants can often appear to be de novo mutations simply because
insufficient sequence coverage was available for one of the parents to detect that the
parent(s) is also a heterozygote (and thus the variant was actually inherited, not
spontaneous).  One simple way to filter such artifacts is to enforce a minimum sequence
depth for each sample.  For example, if we require that at least 50 sequence alignments
were present for mom, dad and child, two of the above variants will be eliminated
as candidates:

.. code-block:: bash

	$ gemini de_novo -d 50 my.db

	family_id	chrom	start	end	ref	alt	gene	impact	impact_severity	in_dbsnp	rs_ids		aaf_1kg_all	aaf_esp_all	clinvar_sig	clinvar_disease_name	clinvar_dbsource	sample1(father)	sample2(mother)	sample3(child; affected)	sample1(depth)	sample2(depth)	sample3(depth)
	1	chr1	17197609	17197610	G	A	BX284668.1	non_syn_coding	MED	1	rs200754171		None	None	None	None	None	G/G	G/G	G/A	104	168	244
	1	chr2	90060872	90060873	A	T	AC009958.1	non_syn_coding	MED	1	rs202041075		None	None	None	None	None	A/A	A/A	A/T	90	238	234
	1	chr3	195505789	195505790	G	C	MUC4	non_syn_coding	MED	1	rs11928301		None	None	None	None	None	G/G	G/G	G/C	250	247	248
	...



============================================================================
``autosomal_recessive``: Find variants meeting an autosomal recessive model.
============================================================================
.. note::

    This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:

    ``gemini load -v my.vcf -p my.ped my.db``

Assuming you have defined the familial relationships between samples when loading
your VCF into GEMINI, one can leverage a built-in tool for identifying variants
that meet an autosomal recessive inheritance pattern. The reported variants
will be restricted to those variants having the potential to impact the
function of affecting protein coding transcripts.

.. code-block:: bash

	$ gemini autosomal_recessive my.db | head

	family_id	chrom	start	end	ref	alt	gene	impact	impact_severity	sample1(father)	sample2(mother)	sample3(child; affected)
	1	chr1	1888192	1888193	C	A	C1orf222	non_syn_coding	MED	C/A	C/A	A/A
	1	chr1	6162053	6162054	T	C	CHD5	non_syn_coding	MED	T/C	T/C	C/C
	1	chr1	6646958	6646968	GCCTGCCTTC	G	ZBTB48	inframe_codon_loss	MED	GCCTGCCTTC/G	GCCTGCCTTC/G	G/G
	1	chr1	11826629	11826630	C	T	C1orf167	non_syn_coding	MED	C/T	C/T	T/T
	1	chr1	11828237	11828238	G	A	C1orf167	non_syn_coding	MED	G/A	G/A	A/A
	1	chr1	11828318	11828319	G	A	C1orf167	non_syn_coding	MED	G/A	G/A	A/A
	1	chr1	11831614	11831615	C	T	C1orf167	non_syn_coding	MED	C/T	C/T	T/T
	1	chr1	11836627	11836628	T	C	C1orf167	non_syn_coding	MED	T/C	T/C	C/C
	1	chr1	11836681	11836682	C	T	C1orf167	non_syn_coding	MED	C/T	C/T	T/T	...



===========================================================================
``autosomal_dominant``: Find variants meeting an autosomal dominant model.
===========================================================================
.. note::

    This tool requires that you identify familial relationships via a PED file
    when loading your VCF into gemini via:

    ``gemini load -v my.vcf -p my.ped my.db``

Assuming you have defined the familial relationships between samples when loading
your VCF into GEMINI, one can leverage a built-in tool for identifying variants
that meet an autosomal dominant inheritance pattern. The reported variants
will be restricted to those variants having the potential to impact the
function of affecting protein coding transcripts.

.. code-block:: bash

	$ gemini autosomal_dominant my.db | head

	family_id	chrom	start	end	ref	alt	gene	impact	impact_severity	sample1(father)	sample2(mother)	sample3(child; affected)
	1	chr1	16855	16856	A	G	WASH7P	splice_donor	HIGH	A/A	A/G	A/G
	1	chr1	881917	881918	G	A	NOC2L	non_syn_coding	MED	G/A	G/G	G/A
	1	chr1	907757	907758	A	G	PLEKHN1	non_syn_coding	MED	A/A	A/G	A/G
	1	chr1	909237	909238	G	C	PLEKHN1	non_syn_coding	MED	G/C	C/C	G/C
	1	chr1	916548	916549	A	G	C1orf170	non_syn_coding	MED	A/G	G/G	A/G
	1	chr1	935221	935222	C	A	HES4	non_syn_coding	MED	C/A	A/A	C/A
	1	chr1	949607	949608	G	A	ISG15	non_syn_coding	MED	G/A	G/G	G/A
	1	chr1	979747	979748	A	T	AGRN	non_syn_coding	MED	A/T	A/A	A/T
	1	chr1	1361529	1361530	C	T	TMEM88B	non_syn_coding	MED	C/T	C/C	C/T



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
