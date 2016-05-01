####################################################################
# 1. Test CADD scores (default)
####################################################################
gemini load -v test.cadd.vcf --skip-gerp-bp test.cadd.db > /dev/null 2>&1

path=`gemini_python ../gemini/anno_info.py`
file="whole_genome_SNVs.tsv.compressed.gz"

if [ -f $path/$file ]; then

echo "    cadd.t01...\c"
echo "chrom	start	end	ref	alt	cadd_raw	cadd_scaled
chr1	30859	30860	G	C	-0.92	0.39
chr1	69269	69270	A	G	-1.34	0.03
chr1	69510	69511	A	G	-1.72	0.01
chr1	874815	874816	C	CT	None	None
chr1	879675	879676	G	A	0.56	7.02
chr1	935491	935492	G	T	0.79	8.17
chr1	1334051	1334057	CTAGAG	C	None	None
chr3	60830533	60830534	A	T	0.18	4.98
chr3	60830762	60830763	G	T	-1.65	0.01" > exp

else

echo "    cadd.t01...\c"
echo "chrom	start	end	ref	alt	cadd_raw	cadd_scaled
chr1	30859	30860	G	C	None	None
chr1	69269	69270	A	G	None	None
chr1	69510	69511	A	G	None	None
chr1	874815	874816	C	CT	None	None
chr1	879675	879676	G	A	None	None
chr1	935491	935492	G	T	None	None
chr1	1334051	1334057	CTAGAG	C	None	None
chr3	60830533	60830534	A	T	None	None
chr3	60830762	60830763	G	T	None	None" > exp

fi

gemini query -q "select chrom, start, end, ref, alt, cadd_raw, cadd_scaled from variants" \
test.cadd.db --header \
> obs
check obs exp
rm obs exp
####################################################################
# 2. Test CADD scores when not loaded with "--skip_cadd"
####################################################################
echo "    cadd.t02...\c"
echo "chrom	start	end	ref	alt	cadd_raw	cadd_scaled
chr1	30859	30860	G	C	None	None
chr1	69269	69270	A	G	None	None
chr1	69510	69511	A	G	None	None
chr1	874815	874816	C	CT	None	None
chr1	879675	879676	G	A	None	None
chr1	935491	935492	G	T	None	None
chr1	1334051	1334057	CTAGAG	C	None	None
chr3	60830533	60830534	A	T	None	None
chr3	60830762	60830763	G	T	None	None" > exp
gemini query -q "select chrom, start, end, ref, alt, cadd_raw, cadd_scaled from variants" \
	   test.no.cadd.db --header \
       > obs
check obs exp
rm obs exp
#####################################################################
