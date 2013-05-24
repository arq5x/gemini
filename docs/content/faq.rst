##############################
F.A.Q.
##############################


=======================================
How can I use PLINK files with GEMINI?
=======================================

Many datasets, especially those derived from GWAS studies, are based on SNP 
genotyping arrays, and are thus stored in the  standard PLINK formats. 
While GEMINI only supports VCF input files, it is relatively straightforward to
convert PLINK datasets to VCF with the PLINK/SEQ toolkit.  We demonstrate a basic
example of the workflow below.

First, load the PLINK BED file into a new PLINK/SEQ project using the instructions
found in the "Load a PLINK binary fileset" section `here <http://atgu.mgh.harvard.edu/plinkseq/input.shtml#plink>`_.

Next, use PLINK/SEQ to convert the project to VCF using the instructions found 
`here <http://atgu.mgh.harvard.edu/plinkseq/output.shtml#vcf>`_.

At this point, you should have a VCF file that is compatible with GEMINI.