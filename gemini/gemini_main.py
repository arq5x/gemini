#!/usr/bin/env python

import os.path
import sys
import argparse
import textwrap
import gemini_load, gemini_query,\
       gemini_region, gemini_stats, gemini_dump, \
       gemini_update

def examples(parser, args):
    
    print
    print "[load] - load a VCF file into a gemini database:"
    print "   gemini load -v my.vcf my.db"
    print "   gemini load -v my.vcf -t snpEff my.db"
    print "   gemini load -v my.vcf -t VEP my.db"
    print

    print "[stats] - report basic statistics about your variants:"
    print "   gemini stats --tstv my.db"
    print "   gemini stats --tstv-coding my.db"
    print "   gemini stats --sfs my.db"
    print "   gemini stats --snp-counts my.db"
    print
    
    print "[query] - explore the database with ad hoc queries:"
    print "   gemini query -q \"select * from variants where is_lof = 1 and aaf <= 0.01\" my.db"
    print "   gemini query -q \"select chrom, pos, gt_bases.NA12878 from variants\" my.db"
    print "   gemini query -q \"select chrom, pos, in_omim, clin_sigs from variants\" my.db"
    print

    print "[dump] - convenient \"data dumps\":"
    print "   gemini dump --variants my.db"
    print "   gemini dump --genotypes my.db"
    print "   gemini dump --samples my.db"
    print
    
    print "[region] - access variants in specific genomic regions:"
    print "   gemini region --reg chr1:100-200 my.db"
    print "   gemini region --gene TP53 my.db"
    print
    
    exit()


def main():
    
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='gemini',
                                     usage='%(prog)s <command> [options]')
    subparsers = parser.add_subparsers(title='commands')

    # $ gemini examples
    parser_examples = subparsers.add_parser('examples', help='show usage examples')
    parser_examples.set_defaults(func=examples)

    # $ gemini load
    parser_load = subparsers.add_parser('load',  help='load a VCF file in gemini database')
    parser_load.add_argument('db', metavar='db',
                            help='The name of the database to be created.')
    parser_load.add_argument('-v', dest='vcf', 
                            help='The VCF file to be loaded.')
    parser_load.add_argument('-t', dest='anno_type', default=None, metavar='STRING',
                             help="The annotations to be used with the input vcf. Options are:\n"
                             "  snpEff  - Annotations as reported by snpEff.\n"
                             "  VEP     - Annotations as reported by VEP.\n"
                             )
    parser_load.add_argument('-p', dest='ped_file', 
                            help='Sample information file in PED+ format.', default=None)
    parser_load.add_argument('--noload-genotypes', dest='noload_genotypes', action='store_true',
                            help='Genotypes exist in the file, but should not be stored.', default=False)
    parser_load.add_argument('--no-genotypes', dest='no_genotypes', action='store_true',
                            help='There are no genotypes in the file (e.g. some 1000G VCFs)', default=False)
    parser_load.set_defaults(func=gemini_load.load)


    # $ gemini query
    parser_query = subparsers.add_parser('query',  help='issue ad hoc SQL queries to the DB')
    parser_query.add_argument('db', metavar='db',  
                            help='The name of the database to be queried.')
    parser_query.add_argument('-q', dest='query', metavar='QUERY_STR',
                            help='The query to be issued to the database')
    parser_query.add_argument('-f', dest='queryfile', metavar='QUERY_FILE', 
                            help='A text file containing the query to be issued.')
    parser_query.add_argument('--header', dest='use_header', action='store_true',
                            help='Add a header of column names to the output.', default=False)
    parser_query.add_argument('--sep', dest='separator', metavar='STRING',
                            help='Output column separator', default="\t")
    parser_query.set_defaults(func=gemini_query.query)


    # $ gemini dump
    parser_dump = subparsers.add_parser('dump', help='shortcuts for extracting data from the DB')
    parser_dump.add_argument('db', metavar='db',  
                            help='The name of the database to be queried.')
    parser_dump.add_argument('--variants', dest='variants', action='store_true',
                            help='Report all rows/columns from the variants table.', default=False)
    parser_dump.add_argument('--genotypes', dest='genotypes', action='store_true',
                            help='Report all rows/columns from the variants table \nwith one line per sample/genotype.', default=False)
    parser_dump.add_argument('--samples', dest='samples', action='store_true',
                            help='Report all rows/columns from the samples table.', default=False)
    parser_dump.add_argument('--header', dest='use_header', action='store_true',
                            help='Add a header of column names to the output.', default=False)
    parser_dump.add_argument('--sep', dest='separator', metavar='STRING',
                            help='Output column separator', default="\t")
    parser_dump.set_defaults(func=gemini_dump.dump)


    # $ gemini region
    parser_region = subparsers.add_parser('region', help='extract variants from specific genomic loci')
    parser_region.add_argument('db', metavar='db',  
                            help='The name of the database to be queried.')
    parser_region.add_argument('--reg', dest='region', metavar='STRING',
                            help='Specify a chromosomal region chr:start-end')
    parser_region.add_argument('--gene', dest='gene', metavar='STRING',
                            help='Specify a gene of interest')
    parser_region.add_argument('--header', dest='use_header', action='store_true',
                            help='Add a header of column names to the output.', default=False)
    parser_region.add_argument('--sep', dest='separator', metavar='STRING',
                            help='Output column separator', default="\t")
    parser_region.set_defaults(func=gemini_region.region)


    # $ gemini stats
    parser_stats = subparsers.add_parser('stats', help='compute useful variant stastics')
    parser_stats.add_argument('db', metavar='db',  
                              help='The name of the database to be queried.')
    parser_stats.add_argument('--tstv', dest='tstv', action='store_true',
                              help='Report the overall ts/tv ratio.', default=False)
    parser_stats.add_argument('--tstv-coding', dest='tstv_coding', action='store_true',
                              help='Report the ts/tv ratio in coding regions.', default=False)
    parser_stats.add_argument('--tstv-noncoding', dest='tstv_noncoding', action='store_true',
                              help='Report the ts/tv ratio in non-coding regions.', default=False)
    parser_stats.add_argument('--snp-counts', dest='snp_counts', action='store_true',
                              help='Report the count of each type of SNP (A->G, G->T, etc.).', default=False)
    parser_stats.add_argument('--sfs', dest='sfs', action='store_true',
                              help='Report the site frequency spectrum of the variants.', default=False)
    parser_stats.add_argument('--mds', dest='mds', action='store_true',
                              help='Report the pairwise genetic distance between the samples.', default=False)
    parser_stats.add_argument('--vars-by-sample', dest='variants_by_sample', action='store_true',
                              help='Report the number of variants observed in each sample.', default=False)
    parser_stats.set_defaults(func=gemini_stats.stats)


    #######################################################
    # "pop update" create the parser for the "update" command
    #######################################################
    parser_get = subparsers.add_parser('update', formatter_class=argparse.RawTextHelpFormatter)
    parser_get.add_argument('db', metavar='db',  help='The name of the database to be created.')
    parser_get.add_argument('-t', dest='table',  help='The table to be updated/modified.')
    parser_get.add_argument('-f', dest='file',   help='The file containing the data to be loaded.')
    parser_get.add_argument('-m', dest='mode',   default="append", 
                           help="Are you appending a new column or updating an existing column?.\n"
                           "Options are:\n"
                           "        append - add a new column to the table (-f).\n"
                           "        update - update an existing column in the table (-f).\n")
    parser_get.set_defaults(func=gemini_update.apply)


    #######################################################
    # parse the args and call the selected function
    #######################################################
    args = parser.parse_args()
    args.func(parser, args)
    
if __name__ == "__main__":
    main()
