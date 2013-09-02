
1. **Return all information from the variant table**

::

  SELECT * from variant


2. **Return all variants with MAF > 0.1**

::

  SELECT * from variant v
  WHERE v.maf > 0.1


3. For each sample, find all genes and variants where the given sample has either:
  a. a single, homozygous, LOF mutation and that sample is the _only_ sample with that mutation.
  b. two or more heterozygous LOF mutations (i.e., a possible compound heterozygote), and that sample is the _only_ sample with those mutations.

::
  SELECT v.*, gt.genotype, g.gene
  FROM variant v,
       genotype gt,
       sample s,
       gene ge
  WHERE v.variant_id = gt.variant_id
    AND gt.sample_id = s.sample_id
    AND v.gene_id    = ge.gene_id
    AND s.sample = "NA00001"
    AND v.lof = 1
    AND (v.num_het = 1 OR v.num_hom_alt = 1)
    AND g.genotype != "hom_ref"
    AND g.genotype != "unknown"