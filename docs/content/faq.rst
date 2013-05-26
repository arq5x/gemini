##############################
F.A.Q.
##############################


========================================
Does GEMINI work with non-human genomes?
========================================

Currently, no.  However, we recognize that the GEMINI framework is suitable to
genetic research in other organisms. This may be a focus of future work.


======================================================
What versions of the human genome does GEMINI support?
======================================================
Currently, we support solely build 37 of the human genome (a.k.a, hg19). We intend to support forthcoming versions of the human genome in future releases.


=======================================
How can I use PLINK files with GEMINI?
=======================================

Many datasets, especially those derived from GWAS studies, are based on SNP 
genotyping arrays, and are thus stored in the  standard PLINK formats. 
While GEMINI only supports VCF input files, it is relatively straightforward to
convert PLINK datasets to VCF with the PLINK/SEQ toolkit.

1. First, load the PLINK BED file into a new PLINK/SEQ project using the instructions
found in the "Load a PLINK binary fileset" section `here <http://atgu.mgh.harvard.edu/plinkseq/input.shtml#plink>`_.

2. Next, use PLINK/SEQ to convert the project to VCF using the instructions found 
`here <http://atgu.mgh.harvard.edu/plinkseq/output.shtml#vcf>`_.

At this point, you should have a VCF file that is compatible with GEMINI.

Alternatively, in his `bcbio <https://github.com/chapmanb/bcbio-nextgen>`_ project, Brad Chapman has written a convenient 
`script <https://github.com/chapmanb/bcbio-nextgen/blob/master/scripts/plink_to_vcf.py>`_ for directly converting PLINK files to VCF.  Below is an example of how to use 
this script.

.. code-block:: bash

	$ plink_to_vcf.py <ped file> <map file> <UCSC reference file in 2bit format)

