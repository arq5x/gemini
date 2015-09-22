check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}
export -f check

<<PEDIGREE
1	1	1_dad	0	0	-1	1
2	1	1_mom	0	0	-1	1
3	1	1_kid	1_dad	1_mom	-1	2
4	2	2_dad	0	0	-1	1
5	2	2_mom	0	0	-1	2
6	2	2_kid	2_dad	2_mom	-1	2
7	3	3_dad	0	0	-1	2
8	3	3_mom	0	0	-1	-9
9	3	3_kid	3_dad	3_mom	-1	2
PEDIGREE

echo "genewise.t1"

echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_filters	n_gene_variants	gene_filters
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	1	2	1
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	1	2	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	1	1" > exp

gemini gene_wise  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
	--gt-filter "(gt_types.1_dad == HOM_REF and gt_types.1_kid == HET and gt_types.1_mom == HOM_REF)" \
    test.auto_dom.db  > obs
check obs exp
rm obs exp

echo "genewise.t2"
echo "ERROR gene-wise: specified --min-filter > the number of --gt-filters" > exp
gemini gene_wise  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
	--gt-filter "(gt_types.1_dad == HOM_REF and gt_types.1_kid == HET and gt_types.1_mom == HOM_REF)" \
	--min-filters 2 \
    test.auto_dom.db  &> obs
check obs exp
rm obs exp

echo "genewise.t3"
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_filters	n_gene_variants	gene_filters
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	1,2	2	1,2
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	1,2	2	1,2
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1,2	1	1,2" > exp
gemini gene_wise  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
	--gt-filter "(gt_types.2_dad == HOM_REF and gt_types.2_kid == HET and gt_types.2_mom == HET)" \
	--gt-filter "(gt_types.1_dad == HOM_REF and gt_types.1_kid == HET and gt_types.1_mom == HOM_REF)" \
	--min-filters 2 \
    test.auto_dom.db  > obs
check obs exp

echo "genewise.t4"
echo "gene	chrom	start	end	ref	alt	impact	impact_severity	variant_filters	n_gene_variants	gene_filters
ASAH2C	chr10	48003991	48003992	C	T	non_syn_coding	MED	1	2	1
ASAH2C	chr10	48004991	48004992	C	T	non_syn_coding	MED	1	2	1
WDR37	chr10	1142207	1142208	T	C	stop_loss	HIGH	1	2	1
WDR37	chr10	1142208	1142209	T	C	stop_loss	HIGH	1	2	1" > exp
gemini gene_wise  \
    --columns "gene, chrom, start, end, ref, alt, impact, impact_severity" \
	--gt-filter  "((gt_types).(phenotype==1).(==HOM_ALT).(none))" \
	--min-filters 1 \
    test.auto_dom.db  > obs
check obs exp
