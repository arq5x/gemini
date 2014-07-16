check()
{
  if diff $1 $2; then
      echo ok
  else
      echo fail
  fi
}
export -f check

########################################################################
# 1. Test default ROH
########################################################################
echo "    roh.t1...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp" > exp
gemini roh test.roh.vcf.db > obs
check obs exp
rm obs exp

########################################################################
# 2. Test ROH with relaxed parameters
########################################################################
echo "    roh.t2...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp
M10475	chr10	1142208	126678092	4	0.0	125535884
M10475	chr10	126678092	172057435	3	0.0001	45379343
M10475	chr11	1142208	126678092	4	0.0	125535884
M10475	chr11	126678092	172057435	3	0.0001	45379343" > exp
gemini roh --min-snps 3 --min-size 10 --min-total-depth 0 test.roh.vcf.db > obs
check obs exp
rm obs exp


########################################################################
# 3. Test ROH with min-depth
########################################################################
echo "    roh.t3...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp
M10475	chr10	1142208	172057435	5	0.0	170915227
M10475	chr11	1142208	172057435	5	0.0	170915227" > exp
gemini roh --min-snps 3 --min-size 10 --min-total-depth 10 test.roh.vcf.db > obs
check obs exp
rm obs exp

########################################################################
# 4. Test ROH with min-depth and lower max-hets
########################################################################
echo "    roh.t4...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp
M10475	chr10	1142208	126678092	3	0.0	125535884
M10475	chr11	1142208	126678092	3	0.0	125535884" > exp
gemini roh --min-snps 3 --min-size 10 --min-total-depth 10 --max-hets 0 test.roh.vcf.db > obs
check obs exp
rm obs exp

########################################################################
# 5. Test ROH raising min depth
########################################################################
echo "    roh.t5...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp
M10475	chr10	1142208	172057435	5	0.0	170915227
M10475	chr11	1142208	172057435	5	0.0	170915227" > exp
gemini roh --min-snps 3 --min-size 10 --min-total-depth 100 --max-hets 0 test.roh.vcf.db > obs
check obs exp
rm obs exp

########################################################################
# 6. Test ROH raising min size 
########################################################################
echo "    roh.t6...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp" > exp
gemini roh --min-snps 3 --min-size 10 --min-total-depth 100 --max-hets 0 --min-size 200000000 test.roh.vcf.db > obs
check obs exp
rm obs exp

########################################################################
# 7. Test ROH raising min size 
########################################################################
echo "    roh.t7...\c"
echo "sample	chrom	run_start	run_end	num_of_snps	density_per_kb	run_length_in_bp
M10475	chr10	1142208	126678092	4	0.0	125535884
M10475	chr11	1142208	126678092	4	0.0	125535884" > exp
gemini roh --min-snps 3 --min-size 10 --min-total-depth 0 --min-size 50000000 test.roh.vcf.db > obs
check obs exp
rm obs exp
