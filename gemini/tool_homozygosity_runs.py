# Report runs of homozygosity for each sample (e.g.hom_ref,hom_alt)
# may consider sample_depth (genotype calls are less likely to be fake),
# allow few hets, unknowns in the region of run
# return runs for a defined window (e.g. >=25 snps for a 500kb region)
import os
import sys
from collections import defaultdict 
from gemini_constants import *
import GeminiQuery

class Site(object):
    def __init__(self, row):
        self.chrom = row['chrom']
        self.end = int(row['end'])
        self.gt_type = None


def sweep_genotypes_for_rohs(args, chrom, samples):
    """
    Sweep through the genotypes for each sample in search of ROHs.

    Note: If the genotype was homozygous, the end position
          of the variant is stored.  Otherwise 'H' for het
          and 'U' for unknown.

    Note: Since we are sweeping through sites until we hit too
    many HETS or UNKNOWNs, we will want to start the next run
    with the last HOM observed in the previous run.
    
    For example (#=hom, H=het, U=unknown).
    if we had the following and min homs is 5 and max hets is 2:
       1 3 7 9 11 H H 17 19 21 23 25
    
    we would report:
      [1 3 7 9 11]
    but we would also want to report:
              [11     17 19 21 23 25]
    """
    hom_count = 0
    het_count = 0
    unk_count = 0
    curr_run = []
    for sample in samples:

        sites = iter(samples[sample])
        for site in sites:

            # retain the last homozygote from previous
            # run. See function docs for details
            if len(curr_run):
                curr_run = [curr_run[-1]]

            # sweep through the active sites until we encounter 
            # too many HETS or UNKNOWN genotypes.
            while het_count <= args.max_hets and unk_count <= args.max_unknowns:
                if site != 'H' and site != 'U':
                    hom_count +=1
                    curr_run.append(site)
                elif site == 'H':
                    het_count += 1
                elif site == 'U':
                    unk_count += 1
                try:
                    site = sites.next()
                except:
                    break

            # skip the current run unless it contains enough sites.
            if len(curr_run) >= args.min_snps:
                run_start = curr_run[0]
                run_end = curr_run[-1]
                run_length = run_end - run_start
                
                # report the run if it is long enough.
                if run_length >= args.min_size:
                    density_per_kb = float (int(len(curr_run)) * 1000) / run_length
                    print "\t".join(str(s) for s in [sample, chrom, 
                        run_start, run_end, 
                        len(curr_run), round(density_per_kb, 2), 
                        run_length])

            # reset for next run
            hom_count = 0
            het_count = 0
            unk_count = 0


def get_homozygosity_runs(args):

    gq = GeminiQuery.GeminiQuery(args.db)
    
    # get a mapping of sample ids to sample indices
    idx2smp = gq.index2sample
    smp2idx = gq.sample2index
    sm_index = []

    # prepare a lookup of just the samples
    # for which the user wishes to search for ROHs
    if args.samples is not None:
        sample_filter = args.samples.strip().split(",")
        for sample in sample_filter:
            sm_index.append(smp2idx[sample])
    else:
        for sample in smp2idx:
            sm_index.append(smp2idx[sample])

    ###########################################################################
    # Phase 1. Retrieve the variants for each chrom/sample
    ###########################################################################
    query  = "SELECT chrom, start, end, gt_types, gt_depths \
              FROM variants \
              WHERE type = 'snp' \
              AND   filter is NULL \
              AND   depth >= " + str(args.min_total_depth) + \
              " ORDER BY chrom, end"

    sys.stderr.write("LOG: Querying and ordering variants by chromosomal position.\n")
    gq.run(query, needs_genotypes=True)

    print "\t".join(['sample','chrom',
        'run_start','run_end',
        'num_of_snps','density_per_kb',
        'run_length_in_bp'])  
    
    variants_seen = 0
    samples = defaultdict(list)
    prev_chrom = None
    curr_chrom = None
    for row in gq:

        variants_seen += 1
        if variants_seen % 10000 == 0:
            sys.stderr.write("LOG: Loaded %d variants. Current variant on %s, position %d.\n" \
                % (variants_seen, row['chrom'], row['end']))

        gt_types = row['gt_types']
        gt_depths = row['gt_depths']
        curr_chrom = row['chrom']

        # the chromosome has changed. search for ROHs in the previous chrom
        if curr_chrom != prev_chrom and prev_chrom is not None:
            sweep_genotypes_for_rohs(args, prev_chrom, samples)
            samples = defaultdict(list)

        # associate the genotype for the variant with each sample
        for idx in sm_index:
            sample = idx2smp[idx]
            gt_type = gt_types[idx]
            depth = gt_depths[idx]

            # the genotype must have had sufficient depth to be considered
            if depth < args.min_genotype_depth:
                continue

            if (gt_type == HOM_ALT or gt_type == HOM_REF):
                samples[sample].append(row['end'])
            elif gt_type == HET:
                samples[sample].append('H')
            elif gt_type == UNKNOWN:
                samples[sample].append('U')

        prev_chrom = curr_chrom

    # search for ROHs in the final chromosome
    sweep_genotypes_for_rohs(args, curr_chrom, samples)


def run(parser, args):

    if os.path.exists(args.db):
        # run the roh caller
        get_homozygosity_runs(args)        

