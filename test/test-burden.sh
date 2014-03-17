###################################################################
# 1. Test high impact burden counts
###################################################################
echo "    burden.t1...\c"
echo "gene	M10475	M10478	M10500	M128215
DHODH	1	0	0	0
WDR37	2	2	2	2
CTBP2	0	0	0	1" > exp
gemini burden test.burden.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test c-alpha test with specified controls and cases
###################################################################
echo "    burden.t2...\c"
echo "gene	T	c	Z	p_value
SYCE1	-0.5	0.25	-1.0	0.841344746069
DHODH	0.0	0.0	nan	nan
WDR37	-1.0	1.5	-0.816496580928	0.792891910879
ASAH2C	-0.5	0.75	-0.57735026919	0.718148569175
CTBP2	0.0	0.0	nan	nan" > exp
gemini burden --controls M10475 M10478 --cases M10500 M128215 --calpha test.burden.db > obs
check obs exp
rm obs exp

####################################################################
## 3. Test c-alpha test with PED file
####################################################################
echo "    burden.t3...\c"
echo "gene	T	c	Z	p_value
SYCE1	-0.5	0.25	-1.0	0.841344746069
DHODH	0.0	0.0	nan	nan
WDR37	-1.0	1.5	-0.816496580928	0.792891910879
ASAH2C	-0.5	0.75	-0.57735026919	0.718148569175
CTBP2	0.0	0.0	nan	nan" > exp
gemini burden --calpha test.burden.db > obs
check obs exp
rm obs exp
