source ./check.sh

gemini query -q "select gt_phred_ll_homref.A, gt_phred_ll_het.A, gt_phred_ll_homalt.A, gt_phred_ll_homref.B, gt_phred_ll_het.B, gt_phred_ll_homalt.B from variants" test.PLs.db > obs
echo "561	0	317	0	99	1176
142	0	251	34	0	311
0	247	2833	744	0	2519
1119	0	3133	438	0	4373
3658	0	2804	1456	0	4749
45	3	0	-1	-1	-1
45	3	0	-1	-1	-1
41	3	0	-1	-1	-1
None	None	None	None	None	None" > exp

check obs exp "GLs.t01"


gemini query -q "select (gt_phred_ll_homref).(*) from variants limit 2" test.PLs.db > obs

echo "561	0	0	0	282	360	0	0	0	137	0	401
142	34	57	404	46	207	99	49	-1	-1	46	37" > exp

check obs exp "GLs.t02"



gemini query --header -q "select (gt_phred_ll_homref).(*) from variants" --gt-filter "(gt_phred_ll_homref).(*).(>=20).(all)" test.PLs.db > obs
echo "gt_phred_ll_homref.A	gt_phred_ll_homref.B	gt_phred_ll_homref.C	gt_phred_ll_homref.D	gt_phred_ll_homref.E	gt_phred_ll_homref.F	gt_phred_ll_homref.G	gt_phred_ll_homref.H	gt_phred_ll_homref.I	gt_phred_ll_homref.J	gt_phred_ll_homref.K	gt_phred_ll_homref.L
1119	438	87	583	383	125	755	88	1649	201	161	103" > exp
check obs exp "GLs.t03"



gemini query --header -q "select (gt_phred_ll_homalt).(*) from variants" --gt-filter "(gt_phred_ll_homalt).(*).(>=20).(all)" test.PLs.db > obs

echo "gt_phred_ll_homalt.A	gt_phred_ll_homalt.B	gt_phred_ll_homalt.C	gt_phred_ll_homalt.D	gt_phred_ll_homalt.E	gt_phred_ll_homalt.F	gt_phred_ll_homalt.G	gt_phred_ll_homalt.H	gt_phred_ll_homalt.I	gt_phred_ll_homalt.J	gt_phred_ll_homalt.K	gt_phred_ll_homalt.L
317	1176	3426	892	254	494	572	491	150	171	748	154
2833	2519	959	2148	1940	1851	192	45	137	1241	1103	840
3133	4373	5113	3872	3687	4458	1869	3878	63	2385	3825	2315
2804	4749	6850	2750	1536	4131	3536	3157	3510	3511	5850	3064" > exp
check obs exp "GLs.t04"


gemini query --header -q "select (gt_phred_ll_het).(*) from variants" --gt-filter "(gt_phred_ll_het).(*).(>=0).(all)" test.PLs.db > obs
echo "gt_phred_ll_het.A	gt_phred_ll_het.B	gt_phred_ll_het.C	gt_phred_ll_het.D	gt_phred_ll_het.E	gt_phred_ll_het.F	gt_phred_ll_het.G	gt_phred_ll_het.H	gt_phred_ll_het.I	gt_phred_ll_het.J	gt_phred_ll_het.K	gt_phred_ll_het.L
0	99	280	78	0	0	48	42	15	0	66	0
247	0	0	184	159	0	0	0	0	105	0	69
0	0	0	0	0	0	0	0	0	0	0	0
0	0	65	0	0	0	0	0	0	0	0	0" > exp
check obs exp "GLs.t05"


gemini query --header -q "select (gt_phred_ll_het).(*) from variants" --gt-filter "(gt_phred_ll_het).(*).(>=0).(any) and gt_types.G == HOM_REF or gt_types.G == HET" test.PLs.db > obs
echo "gt_phred_ll_het.A	gt_phred_ll_het.B	gt_phred_ll_het.C	gt_phred_ll_het.D	gt_phred_ll_het.E	gt_phred_ll_het.F	gt_phred_ll_het.G	gt_phred_ll_het.H	gt_phred_ll_het.I	gt_phred_ll_het.J	gt_phred_ll_het.K	gt_phred_ll_het.L
0	99	280	78	0	0	48	42	15	0	66	0
0	0	0	0	0	0	0	0	-1	-1	0	0
247	0	0	184	159	0	0	0	0	105	0	69
0	0	0	0	0	0	0	0	0	0	0	0
0	0	65	0	0	0	0	0	0	0	0	0" > exp
check obs exp "GLs.t06"
