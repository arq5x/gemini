#!/usr/bin/env python

import os.path
import sys
import argparse
import textwrap
import gemini_load, gemini_get, gemini_update

def main():
    
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='gemini')
    subparsers = parser.add_subparsers()


    #######################################################
    # "pop load": create the parser for the "load" command
    #######################################################
    parser_load = subparsers.add_parser('load')
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


    #######################################################
    # "pop get" create the parser for the "get" command
    #######################################################
    parser_get = subparsers.add_parser('get', formatter_class=argparse.RawTextHelpFormatter)
    parser_get.add_argument('db', metavar='db',  
                            help='The name of the database to be queried.')
    parser_get.add_argument('-q', dest='query', metavar='STRING',
                            help='The query to be issued to the database')
    parser_get.add_argument('-f', dest='queryfile', metavar='FILE', 
                            help='A text file containing the query to be issued.')
    parser_get.add_argument('--header', dest='use_header', action='store_true',
                            help='Add a header of column names to the output.', default=False)
    parser_get.add_argument('--sep', dest='separator', metavar='STRING',
                            help='Output column separator', default="\t")
    parser_get.add_argument('--reg', dest='region', metavar='STRING',
                            help='Specify a chromosomal region chr:start-end')
    parser_get.add_argument('--gene', dest='gene', metavar='STRING',
                            help='Specify a gene of interest')
    parser_get.add_argument('--prec', dest='precision', metavar='INT',
                            help='Output precision. Only relevant to certain shortcuts (e.g. sfs)', default=3)
    parser_get.add_argument('-s', dest='shortcut', default="variants", metavar='STRING',
                            help="The shortcut query to be issued to the database. Options are:\n"
                           "  variants        - all rows/columns from the variants table.\n"
                           "  samples         - all rows/columns from the samples table.\n"
                           "  genotypes       - one line per variant and sample genotype.\n"
                           "  region          - returns all variants in a given region (--reg chr:start-end).\n"
                           "  gene            - returns all variants affecting a specific gene (e.g., TP53).\n"
                           "  tstv            - report the transition/transversion (Ts/Tv) ratio.\n"
                           "  tstv-coding     - report the Ts/Tv ratio in coding regions.\n"
                           "  tstv-noncoding  - report the Ts/Tv ratio in non-coding regions.\n"
                           "  snp-counts      - the count of each type of SNP (A->G, G->T, etc.).\n"
                           "  sfs             - report the site frequency spectrum of the variants.\n"
                           "  mds             - returns the genetic distance between any two samples.\n"
                           )
    parser_get.set_defaults(func=gemini_get.get)


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
