#############################################################################
# 1. Test pfam domain track in variants table                                         
#############################################################################
echo "   domain.t01...\c"                                                   
echo "chr1	30859	30860	None
chr1	69269	69270	7tm_1
chr1	69510	69511	7tm_1
chr1	874815	874816	None
chr1	879675	879676	None
chr1	935491	935492	None
chr1	1334051	1334057	None
chr1	31896607	31896608	Serinc" > exp                                                                              
gemini query -q "select chrom, start, end, pfam_domain from variants" test3.snpeff.db > obs
check obs exp                                                               
rm obs exp                                                                   
