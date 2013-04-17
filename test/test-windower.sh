#############################################################################
# 1. Test window data for hwe (mean)                                         
#############################################################################
echo "    window.t01...\c"                                                   
echo "chr1	30000	40000	0.856
chr1	60000	70000	0.7865" > exp                                                                              
gemini windower -w 10000 -t hwe -o mean \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 2. Test window data for hwe (median)                                       
#############################################################################
# echo "    window.t02...\c"                                                   
# echo "chr1	30000	40000	0.7751
# chr1	60000	70000	0.1573" > exp                                        
#                                                                              
# gemini windower -w 10000 -t hwe -o median \
#               test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
# check obs exp                                                                
# rm obs exp                                                                   
                                                                             
#############################################################################
# 3. Test window data for hwe (min)                                          
#############################################################################
echo "    window.t03...\c"                                                   
echo "chr1	30000	40000	0.505
chr1	60000	70000	0.1573" > exp                                        
                                                                             
gemini windower -w 10000 -t hwe -o min \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 4. Test window data for hwe (max)                                          
#############################################################################
echo "    window.t04...\c"                                                   
echo "chr1	30000	40000	1.0
chr1	60000	70000	1.0" > exp                                           
                                                                             
gemini windower -w 10000 -t hwe -o max \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 5. Test window data for hwe (collapse)                                     
#############################################################################
echo "    window.t05...\c"                                                   
echo "chr1	30000	40000	1.000000,1.000000,0.775097,0.504985,1.000000
chr1	60000	70000	1.000000,1.000000,0.157299,0.775097,1.000000" > exp  
                                                                             
gemini windower -w 10000 -t hwe -o collapse \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 6. Test window data for nucl_div (mean)                                    
#############################################################################
echo "    window.t06...\c"                                                   
echo "chr1	30000	31000	0.1357
chr1	69000	70000	0.1833" > exp                                        
                                                                             
gemini windower -w 1000 -t nucl_div -o mean \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp
rm obs exp                                                                   
                                                                             
#############################################################################
# 7. Test window data for nucl_diversity (median)                            
#############################################################################
# echo "    window.t07...\c"                                                   
# echo "chr1	30000	31000	0.25
# chr1	69000	70000	0.6667" > exp                                        
#                                                                              
# gemini windower -w 1000 -t nucl_div -o median \
#               test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
# check obs exp                                                                
# rm obs exp                                                                   
                                                                             
#############################################################################
# 8. Test window data for nucl_diversity (min)                               
#############################################################################
echo "    window.t08...\c"                                                   
echo "chr1	30000	31000	0.0
chr1	69000	70000	0.0" > exp                                           
                                                                             
gemini windower -w 1000 -t nucl_div -o min \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 9. Test window data for nucl_div (max)                                     
#############################################################################
echo "    window.t09...\c"                                                   
echo "chr1	30000	31000	0.4286
chr1	69000	70000	0.6667" > exp                                        
                                                                             
gemini windower -w 1000 -t nucl_div -o max \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 10. Test window data for nucl_diversity (collapse)                         
#############################################################################
echo "    window.t10...\c"                                                   
echo "chr1	30000	31000	0.000000,0.000000,0.250000,0.428571,0.000000
chr1	69000	70000	0.000000,0.000000,0.666667,0.250000,0.000000" > exp  
                                                                             
gemini windower -w 1000 -t nucl_div -o collapse \
              test.snpeff.vcf.db | awk '{if($4!=".")print}' > obs            
check obs exp                                                                
rm obs exp                                                                   
                                                                             
#############################################################################
# 11. Test window data for -s option (default -t(hwe), -o(mean))             
#############################################################################
# echo "    window.t11...\c"
# echo "chr1	30000	31000   0.856
# chr1	30500	31500   0.856
# chr1    69000   70000   0.7865
# chr1    69500   70500   0.6441" > exp
# 
# gemini windower -w 1000 -s 500 \
#               test.snpeff.vcf.db | awk '{if($4!=".")print}' \
#     > obs
# check obs exp
# rm obs exp
