############################
Using the GEMINI API
############################

.. automodule:: gemini

The GeminiQuery class
=====================
.. autoclass:: GeminiQuery
   :members: run, header, sample2index, index2sample
   :undoc-members:




============================================
Extracting the VCF INFO tags with GEMINI API
============================================
The GEMINI API is useful to extract the individual tags within the INFO field of a VCF (stored as a 
compressed dictionary in the variants table). This would be of particular interest to those who want 
to add custom annotations to their VCF and still be able to access the individual tags programmatically. 
Here is an example where we try to extract the dbNSFP fields from the 'INFO' tag of a VCF, using the API.

.. code-block:: bash
	
	#!/usr/bin/env python
	import sys
	from gemini import GeminiQuery
 
	database = sys.argv[1]
 	gq = GeminiQuery(database)
	query = "SELECT variant_id, chrom, start, end, ref, alt, info \
	         FROM variants"
 
	gq.run(query)
 
	for row in gq:
	    try:
		print "\t".join([str(row['chrom']), str(row['start']), str(row['end']),
			      str(row['ref']), str(row['alt']), str(row.info['dbNSFP_SIFT_pred'])])
	    except KeyError:
	        pass
			
	# yields		
	chr1	906272	906273	C	T	P|D|P
	chr1	906273	906274	C	A	D|D|D
	chr1	906276	906277	T	C	D|D|D
	chr1	906297	906298	G	T	B|B|B
	chr1	1959074	1959075	A	C	D
	chr1	1959698	1959699	G	A	B
	chr1	1961452	1961453	C	T	P
	chr1	2337953	2337954	C	T	D