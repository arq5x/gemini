#!/usr/bin/env python
import os.path
import sys
import argparse
import gemini.version


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

    print "[tools] - there are also many specific tools available"
    print "   1. Find compound heterozygotes."
    print "     gemini comp_hets my.db"
    print

    exit()

def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='gemini', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed gemini version",
                        action="version",
                        version="%(prog)s " + str(gemini.version.__version__))
    parser.add_argument('--annotation-dir', dest='annotation_dir',
                             help='Path to the annotation database.\n'
                                'This argument is optional and if given will take precedence over the default location stored in the gemini config file.')
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command')

    #########################################
    # $ gemini examples
    #########################################
    parser_examples = subparsers.add_parser('examples',
                                            help='show usage examples')
    parser_examples.set_defaults(func=examples)

    #########################################
    # $ gemini load
    #########################################
    parser_load = subparsers.add_parser('load',
                                        help='load a VCF file in gemini database')
    parser_load.add_argument('db', metavar='db',
                             help='The name of the database to be created.')
    parser_load.add_argument('-v', dest='vcf',
                             help='The VCF file to be loaded.')
    parser_load.add_argument('-t', dest='anno_type',
                             default=None, choices=["snpEff", "VEP"],
                             help="The annotations to be used with the input vcf.")
    parser_load.add_argument('-p', dest='ped_file',
                             help='Sample information file in PED+ format.',
                             default=None)
    parser_load.add_argument('--skip-gerp-bp',
                             dest='skip_gerp_bp',
                             action='store_true',
                             help='Do not load GERP scores at base pair resolution. Loaded by default.',
                             default=False)
    parser_load.add_argument('--skip-cadd',
                             dest='skip_cadd',
                             action='store_true',
                             help='Do not load CADD scores. Loaded by default',
                             default=False)
    parser_load.add_argument('--skip-gene-tables',
                             dest='skip_gene_tables',
                             action='store_true',
                             help='Do not load gene tables. Loaded by default.',
                             default=False)
    parser_load.add_argument('--skip-info-string',
                             dest='skip_info_string',
                             action='store_true',
                             help='Do not load INFO string from VCF file to reduce DB size. Loaded by default',
                             default=False)
    parser_load.add_argument('--no-load-genotypes',
                             dest='no_load_genotypes',
                             action='store_true',
                             help='Genotypes exist in the file, but should not be stored.',
                             default=False)
    parser_load.add_argument('--no-genotypes',
                             dest='no_genotypes',
                             action='store_true',
                             help='There are no genotypes in the file (e.g. some 1000G VCFs)',
                             default=False)
    parser_load.add_argument('--cores', dest='cores',
                             default=1,
                             type=int,
                             help="Number of cores to use to load in parallel.")
    parser_load.add_argument('--scheduler', dest='scheduler', default=None,
                             choices=["lsf", "sge", "slurm", "torque"],
                             action=IPythonAction,
                             help='Cluster scheduler to use.')
    parser_load.add_argument('--queue', dest='queue',
                             default=None, help='Cluster queue to use.')
    parser_load.add_argument('--passonly',
                             dest='passonly',
                             default=False,
                             action='store_true',
                             help="Keep only variants that pass all filters.")
    parser_load.add_argument('--test-mode',
                         dest='test_mode',
                         action='store_true',
                         help='Load in test mode (faster)',
                         default=False)

    def load_fn(parser, args):
        import gemini_load
        gemini_load.load(parser, args)

    parser_load.set_defaults(func=load_fn)
    #########################################
    # $ gemini amend
    #########################################
    parser_amend = subparsers.add_parser('amend',
                                         help="Amend an already loaded GEMINI database.")
    parser_amend.add_argument('db',
                              metavar='db',
                              help='The name of the database to be amended.')
    parser_amend.add_argument('--sample',
                              metavar='sample',
                              default=None,
                              help='New sample information file to load')
    def amend_fn(parser, args):
        import gemini_amend
        gemini_amend.amend(parser, args)
    parser_amend.set_defaults(func=amend_fn)

    #########################################
    # $ gemini load_chunk
    #########################################
    parser_loadchunk = subparsers.add_parser('load_chunk',
                                             help='load a VCF file in gemini database')
    parser_loadchunk.add_argument('db',
                                  metavar='db',
                                  help='The name of the database to be created.')
    parser_loadchunk.add_argument('-v',
                                  dest='vcf',
                                  help='The VCF file to be loaded.')
    parser_loadchunk.add_argument('-t',
                                  dest='anno_type',
                                  default=None,
                                  metavar='STRING',
                                  help="The annotations to be used with the input vcf. Options are:\n"
                                  "  snpEff  - Annotations as reported by snpEff.\n"
                                  "  VEP     - Annotations as reported by VEP.\n"
                                  )
    parser_loadchunk.add_argument('-o',
                                  dest='offset',
                                  help='The starting number for the variant_ids',
                                  default=None)
    parser_loadchunk.add_argument('-p',
                                  dest='ped_file',
                                  help='Sample information file in PED+ format.',
                                  default=None)
    parser_loadchunk.add_argument('--no-load-genotypes',
                                  dest='no_load_genotypes',
                                  action='store_true',
                                  help='Genotypes exist in the file, but should not be stored.',
                                  default=False)
    parser_loadchunk.add_argument('--no-genotypes',
                                  dest='no_genotypes',
                                  action='store_true',
                                  help='There are no genotypes in the file (e.g. some 1000G VCFs)',
                                  default=False)
    parser_loadchunk.add_argument('--skip-gerp-bp',
                                  dest='skip_gerp_bp',
                                  action='store_true',
                                  help='Do not load GERP scores at base pair resolution. Loaded by default.',
                                  default=False)
    parser_loadchunk.add_argument('--skip-cadd',
                                 dest='skip_cadd',
                                 action='store_true',
                                 help='Do not load CADD scores. Loaded by default',
                                 default=False)
    parser_loadchunk.add_argument('--skip-gene-tables',
                             dest='skip_gene_tables',
                             action='store_true',
                             help='Do not load gene tables. Loaded by default.',
                             default=False)
    parser_loadchunk.add_argument('--skip-info-string',
                                  dest='skip_info_string',
                                  action='store_true',
                                  help='Do not load INFO string from VCF file to reduce DB size. Loaded by default',
                                  default=False)
    parser_loadchunk.add_argument('--passonly',
                                  dest='passonly',
                                  default=False,
                                  action='store_true',
                                  help="Keep only variants that pass all filters.")
    parser_loadchunk.add_argument('--test-mode',
                         dest='test_mode',
                         action='store_true',
                         help='Load in test mode (faster)',
                         default=False)
    def loadchunk_fn(parser, args):
        import gemini_load_chunk
        gemini_load_chunk.load(parser, args)
    parser_loadchunk.set_defaults(func=loadchunk_fn)

    #########################################
    # $ gemini merge_chunks
    #########################################
    parser_mergechunks = subparsers.add_parser('merge_chunks',
            help='combine intermediate db files into the final gemini ')
    parser_mergechunks.add_argument('--db',
            dest='db',
            help='The name of the final database to be loaded.')
    parser_mergechunks.add_argument('--chunkdb',
            nargs='*',
            dest='chunkdbs',
            action='append')

    def mergechunk_fn(parser, args):
        import gemini_merge_chunks
        gemini_merge_chunks.merge_chunks(parser, args)
    parser_mergechunks.set_defaults(func=mergechunk_fn)

    #########################################
    # $ gemini query
    #########################################
    parser_query = subparsers.add_parser('query',
            help='issue ad hoc SQL queries to the DB')
    parser_query.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_query.add_argument('-q',
            dest='query',
            metavar='QUERY_STR',
            help='The query to be issued to the database')
    parser_query.add_argument('--gt-filter',
            dest='gt_filter',
            metavar='STRING',
            help='Restrictions to apply to genotype values')
    parser_query.add_argument('--show-samples',
                              dest='show_variant_samples',
                              action='store_true',
                              default=False,
                              help=('Add a column of all sample names with a variant to each '
                                    'variant.'))
    parser_query.add_argument('--show-families',
                              dest='show_families',
                              action='store_true',
                              default=False,
                              help=('Add a column listing all of the families '
                                    'with a variant to each variant.'))
    parser_query.add_argument('--family-wise',
                              dest='family_wise',
                              default=False,
                              action='store_true',
                              help=('Perform the sample-filter on a family-wise '
                                    'basis.'))
    parser_query.add_argument('--min-kindreds',
                              dest='min_kindreds',
                              default=1,
                              type=int,
                              help=('Minimum number of families for a variant passing '
                                    'a family-wise filter to be in.'))
    parser_query.add_argument('--sample-delim',
            dest='sample_delim',
            metavar='STRING',
            help='The delimiter to be used with the --show-samples option.',
            default=',')

    parser_query.add_argument('--header',
            dest='use_header',
            action='store_true',
            help='Add a header of column names to the output.',
            default=False)
    parser_query.add_argument('--sample-filter',
                              dest='sample_filter',
                              help='SQL filter to use to filter the sample table',
                              default=None)
    parser_query.add_argument('--in',
                              dest='in_subject',
                              nargs='*',
                              help=('A variant must be in either all, none or any '
                                    'samples passing the --sample-query filter.'),
                              choices=['all', 'none', 'any', 'only'],
                              default=['any'])
    parser_query.add_argument('--format',
                              dest='format',
                              default='default',
                              help='Format of output (JSON, TPED or default)')
    parser_query.add_argument('--region',
                              dest='region',
                              default=None,
                              help=('Restrict query to this region, '
                                    'e.g. chr1:10-20.'))
    parser_query.add_argument('--carrier-summary-by-phenotype',
                              dest='carrier_summary',
                              default=None,
                              help=('Output columns of counts of carriers and '
                                    'non-carriers stratified by the given '
                                    'sample phenotype column'))
    parser_query.add_argument('--dgidb',
                              dest='dgidb',
                              action='store_true',
                              help='Request drug-gene interaction info from DGIdb.',
                              default=False)
    def query_fn(parser, args):
        import gemini_query
        gemini_query.query(parser, args)

    parser_query.set_defaults(func=query_fn)

    #########################################
    # $ gemini dump
    #########################################
    parser_dump = subparsers.add_parser('dump',
            help='shortcuts for extracting data from the DB')
    parser_dump.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_dump.add_argument('--variants',
            dest='variants',
            action='store_true',
            help='Report all rows/columns from the variants table.',
            default=False)
    parser_dump.add_argument('--genotypes',
            dest='genotypes',
            action='store_true',
            help='Report all rows/columns from the variants table \nwith one line per sample/genotype.',
            default=False)
    parser_dump.add_argument('--samples',
            dest='samples',
            action='store_true',
            help='Report all rows/columns from the samples table.',
            default=False)
    parser_dump.add_argument('--header',
            dest='use_header',
            action='store_true',
            help='Add a header of column names to the output.',
            default=False)
    parser_dump.add_argument('--sep',
            dest='separator',
            metavar='STRING',
            help='Output column separator',
            default="\t")
    parser_dump.add_argument('--tfam',
                             dest='tfam',
                             action='store_true',
                             default=False,
                             help='Output sample information to TFAM format.')
    def dump_fn(parser, args):
        import gemini_dump
        gemini_dump.dump(parser, args)
    parser_dump.set_defaults(func=dump_fn)

    #########################################
    # $ gemini region
    #########################################
    parser_region = subparsers.add_parser('region',
            help='extract variants from specific genomic loci')
    parser_region.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_region.add_argument('--reg',
            dest='region',
            metavar='STRING',
            help='Specify a chromosomal region chr:start-end')
    parser_region.add_argument('--gene',
            dest='gene',
            metavar='STRING',
            help='Specify a gene of interest')
    parser_region.add_argument('--header',
            dest='use_header',
            action='store_true',
            help='Add a header of column names to the output.',
            default=False)
    parser_region.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_region.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_region.add_argument('--show-samples',
                               dest='show_variant_samples',
                               action='store_true',
                               default=False,
                                help=('Add a column of all sample names with a variant to each '
                                      'variant.'))
    parser_region.add_argument('--format',
                              dest='format',
                              default='default',
                              help='Format of output (JSON, TPED or default)')
    def region_fn(parser, args):
        import gemini_region
        gemini_region.region(parser, args)
    parser_region.set_defaults(func=region_fn)

    #########################################
    # $ gemini stats
    #########################################
    parser_stats = subparsers.add_parser('stats',
            help='compute useful variant stastics')
    parser_stats.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_stats.add_argument('--tstv',
            dest='tstv',
            action='store_true',
            help='Report the overall ts/tv ratio.',
            default=False)
    parser_stats.add_argument('--tstv-coding',
            dest='tstv_coding',
            action='store_true',
            help='Report the ts/tv ratio in coding regions.',
            default=False)
    parser_stats.add_argument('--tstv-noncoding',
            dest='tstv_noncoding',
            action='store_true',
            help='Report the ts/tv ratio in non-coding regions.',
            default=False)
    parser_stats.add_argument('--snp-counts',
            dest='snp_counts',
            action='store_true',
            help='Report the count of each type of SNP (A->G, G->T, etc.).',
            default=False)
    parser_stats.add_argument('--sfs',
            dest='sfs',
            action='store_true',
            help='Report the site frequency spectrum of the variants.',
            default=False)
    parser_stats.add_argument('--mds',
            dest='mds',
            action='store_true',
            help='Report the pairwise genetic distance between the samples.',
            default=False)
    parser_stats.add_argument('--vars-by-sample',
            dest='variants_by_sample',
            action='store_true',
            help='Report the number of variants observed in each sample.',
            default=False)
    parser_stats.add_argument('--gts-by-sample',
            dest='genotypes_by_sample',
            action='store_true',
            help='Report the count of each genotype class obs. per sample.',
            default=False)
    parser_stats.add_argument('--summarize',
            dest='query',
            metavar='QUERY_STR',
            default=None,
            help='The query to be issued to the database to summarize')
    parser_stats.add_argument('--gt-filter',
            dest='gt_filter',
            metavar='STRING',
            help='Restrictions to apply to genotype values')
    def stats_fn(parser, args):
        import gemini_stats
        gemini_stats.stats(parser, args)
    parser_stats.set_defaults(func=stats_fn)

    #########################################
    # gemini annotate
    #########################################
    parser_get = subparsers.add_parser('annotate',
            help='Add new columns for custom annotations')
    parser_get.add_argument('db',
            metavar='db',
            help='The name of the database to be updated.')
    parser_get.add_argument('-f',
            dest='anno_file',
            help='The TABIX\'ed BED file containing the annotations')
    parser_get.add_argument('-c',
            dest='col_names',
            help='The name(s) of the column(s) to be added to the variant table.')
    parser_get.add_argument('-a',
            dest='anno_type',
            help='How should the annotation file be used? (def. extract)',
            default="extract",
            choices=['boolean', 'count', 'extract'])
    parser_get.add_argument('-e',
            dest='col_extracts',
            help='Column(s) to extract information from for list annotations.')
    parser_get.add_argument('-t',
            dest='col_types',
            help='What data type(s) should be used to represent the new values '
                 'in the database? '
                 'Any of {integer, float, text}')
    parser_get.add_argument('-o',
            dest='col_operations',
            help='Operation(s) to apply to the extract column values '
                  'in the event that a variant overlaps multiple annotations '
                  'in your annotation file (-f).'
                  'Any of {mean, median, min, max, mode, list, uniq_list, first, last}')
    def annotate_fn(parser, args):
        import gemini_annotate
        gemini_annotate.annotate(parser, args)
    parser_get.set_defaults(func=annotate_fn)

    #########################################
    # gemini windower
    #########################################
    parser_get = subparsers.add_parser('windower',
            help='Compute statistics across genome \"windows\"')
    parser_get.add_argument('db',
            metavar='db',
            help='The name of the database to be updated.')
    parser_get.add_argument('-w',
            dest='window_size',
            default=1000000,
            help='The name of the column to be added to the variant table.')
    parser_get.add_argument('-s',
            dest='step_size',
            default=0,
            help="The step size for the windows in bp.\n")
    parser_get.add_argument('-t',
            dest='analysis_type',
            help='The type of windowed analysis requested.',
            choices=['nucl_div', 'hwe'],
            default='hwe')
    parser_get.add_argument('-o',
            dest='op_type',
            help='The operation that should be applied to the -t values.',
            choices=['mean', 'median', 'min', 'max', 'collapse'],
            default='mean')
    def windower_fn(parser, args):
        import gemini_windower
        gemini_windower.windower(parser, args)
    parser_get.set_defaults(func=windower_fn)

    #########################################
    # gemini db_info
    #########################################
    parser_get = subparsers.add_parser('db_info',
            help='Get the names and types of cols. database tables')
    parser_get.add_argument('db',
            metavar='db',
            help='The name of the database to be updated.')
    def dbinfo_fn(parser, args):
        import gemini_dbinfo
        gemini_dbinfo.db_info(parser, args)
    parser_get.set_defaults(func=dbinfo_fn)

    #########################################
    # $ gemini comp_hets
    #########################################
    parser_comp_hets = subparsers.add_parser('comp_hets',
            help='Identify compound heterozygotes')
    parser_comp_hets.add_argument('db',
            metavar='db',
            help='The name of the database to be created.')
    parser_comp_hets.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_comp_hets.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_comp_hets.add_argument('--only-affected',
            dest='only_affected',
            action='store_true',
            help='Report solely those compund heterozygotes impacted a sample \
                  labeled as affected.',
            default=False)
    parser_comp_hets.add_argument('--ignore-phasing',
            dest='ignore_phasing',
            action='store_true',
            help='Ignore phasing when screening for compound hets. \
                  Candidates are inherently _putative_.',
            default=False)
    def comp_hets_fn(parser, args):
        import tool_compound_hets
        tool_compound_hets.run(parser, args)
    parser_comp_hets.set_defaults(func=comp_hets_fn)

    #########################################
    # $ gemini pathways
    #########################################
    parser_pathway = subparsers.add_parser('pathways',
            help='Map genes and variants to KEGG pathways')
    parser_pathway.add_argument('db',
            metavar='db',
            help='The name of the database to be queried')
    parser_pathway.add_argument('-v',
            dest='version',
            default='68',
            metavar='STRING',
            help="Version of ensembl genes to use. "
                 "Supported versions: 66 to 71\n"
            )
    parser_pathway.add_argument('--lof',
            dest='lof',
            action='store_true',
            help='Report pathways for indivs/genes/sites with LoF variants',
            default=False)
    def pathway_fn(parser, args):
        import tool_pathways
        tool_pathways.pathways(parser, args)
    parser_pathway.set_defaults(func=pathway_fn)

    #########################################
    # $ gemini lof_sieve
    #########################################
    parser_lof_sieve = subparsers.add_parser('lof_sieve',
            help='Prioritize LoF mutations')
    parser_lof_sieve.add_argument('db',
            metavar='db',
            help='The name of the database to be queried')
    def lof_sieve_fn(parser, args):
        import tool_lof_sieve
        tool_lof_sieve.lof_sieve(parser, args)
    parser_lof_sieve.set_defaults(func=lof_sieve_fn)

    #########################################
    # $ gemini burden
    #########################################
    burden_help = ("Gene-level genetic burden tests. By default counts all "
                   "variants with high impact in coding regions "
                   "as contributing to burden.")

    parser_burden = subparsers.add_parser('burden',
                                          help=burden_help)
    parser_burden.add_argument('--nonsynonymous', action='store_true',
                               default=False,
                               help=("Count all nonsynonymous variants as "
                                     "contributing burden."))
    parser_burden.add_argument('--cases',
                               dest='cases',
                               nargs='*',
                               help=('Space separated list of cases for '
                                     'association testing.'))
    parser_burden.add_argument('--controls',
                               nargs='*',
                               dest='controls',
                               help=('Space separated list of controls for '
                                     'association testing.'))
    parser_burden.add_argument('--calpha',
                               action='store_true',
                               default=False,
                               help="Run the C-alpha association test.")
    parser_burden.add_argument('--permutations',
                               default=0,
                               type=int,
                               help=("Number of permutations to run for the "
                                     "C-alpha test (try 1000 to start)."))
    parser_burden.add_argument('--min-aaf',
                               dest='min_aaf',
                               type=float,
                               default=0.0,
                               help='The min. alt. allele frequency for a '
                                     'variant to be included.')
    parser_burden.add_argument('--max-aaf',
                               dest='max_aaf',
                               type=float,
                               default=1.0,
                               help='The max. alt. allele frequency for a '
                                     'variant to be included.')
    parser_burden.add_argument('--save_tscores', default=False,
                               action='store_true',
                               help='Save the permuted T-scores to a file.')
    parser_burden.add_argument('db',
                               metavar='db',
                               help='The name of the database to be queried.')

    def burden_fn(parser, args):
        import tool_burden_tests
        tool_burden_tests.burden(parser, args)
    parser_burden.set_defaults(func=burden_fn)

    #########################################
    # $ gemini interactions
    #########################################
    parser_interaction = subparsers.add_parser('interactions',
            help='Find interaction partners for a gene in sample variants(default mode)')
    parser_interaction.add_argument('db',
            metavar='db',
            help='The name of the database to be queried')
    parser_interaction.add_argument('-g',
            dest='gene',
            help='Gene to be used as a root in BFS/shortest_path')
    parser_interaction.add_argument('-r',
            dest='radius',
            type=int,
            help="Set filter for BFS:\n"
                 "valid numbers starting from 0")
    parser_interaction.add_argument('--var',
            dest='var_mode',
            help='var mode: Returns variant info (e.g. impact, biotype) for interacting genes',
            action='store_true',
            default=False)
    def interactions_fn(parser, args):
        import tool_interactions
        tool_interactions.genequery(parser, args)
    parser_interaction.set_defaults(func=interactions_fn)

    #########################################
    # gemini lof_interactions
    #########################################
    parser_interaction = subparsers.add_parser('lof_interactions',
            help='Find interaction partners for a lof gene in sample variants(default mode)')
    parser_interaction.add_argument('db',
            metavar='db',
            help='The name of the database to be queried')
    parser_interaction.add_argument('-r',
            dest='radius',
            type=int,
            help="set filter for BFS:\n")
    parser_interaction.add_argument('--var',
            dest='var_mode',
            help='var mode: Returns variant info (e.g. impact, biotype) for interacting genes',
            action='store_true',
            default=False)
    def lof_interactions_fn(parser, args):
        import tool_interactions
        tool_interactions.lofgenequery(parser, args)
    parser_interaction.set_defaults(func=lof_interactions_fn)

    #########################################
    # $ gemini autosomal_recessive
    #########################################
    parser_auto_rec = subparsers.add_parser('autosomal_recessive',
            help='Identify variants meeting an autosomal \
                  recessive inheritance model')
    parser_auto_rec.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_auto_rec.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_auto_rec.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_auto_rec.add_argument('--min-kindreds',
            dest='min_kindreds',
            type=int,
            default=1,
            help='The min. number of kindreds that must have a candidate variant in a gene.')
    parser_auto_rec.add_argument('-d',
            dest='min_sample_depth',
            type=int,
            help="The minimum aligned\
              sequence depth (genotype DP) req'd for\
              each sample (def. = 0)",
            default=0)
    def autosomal_recessive_fn(parser, args):
        import tool_autosomal_recessive
        tool_autosomal_recessive.run(parser, args)
    parser_auto_rec.set_defaults(func=autosomal_recessive_fn)

    #########################################
    # $ gemini autosomal_dominant
    #########################################
    parser_auto_dom = subparsers.add_parser('autosomal_dominant',
            help='Identify variants meeting an autosomal \
                  dominant inheritance model')
    parser_auto_dom.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_auto_dom.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_auto_dom.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_auto_dom.add_argument('--min-kindreds',
            dest='min_kindreds',
            type=int,
            default=1,
            help='The min. number of kindreds that must have a candidate variant in a gene.')
    parser_auto_dom.add_argument('-d',
            dest='min_sample_depth',
            type=int,
            help="The minimum aligned\
              sequence depth (genotype DP) req'd for\
              each sample (def. = 0)",
            default=0)
    def autosomal_dominant_fn(parser, args):
        import tool_autosomal_dominant
        tool_autosomal_dominant.run(parser, args)
    parser_auto_dom.set_defaults(func=autosomal_dominant_fn)

    #########################################
    # $ gemini de_novo
    #########################################
    parser_de_novo = subparsers.add_parser('de_novo',
            help='Identify candidate de novo mutations')
    parser_de_novo.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_de_novo.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_de_novo.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_de_novo.add_argument('-d',
            dest='min_sample_depth',
            type=int,
            help="The minimum aligned\
                  sequence depth (genotype DP) req'd for\
                  each sample (def. = 0)",
            default=0)
    def de_novo_fn(parser, args):
        import tool_de_novo_mutations
        tool_de_novo_mutations.run(parser, args)
    parser_de_novo.set_defaults(func=de_novo_fn)


    #########################################
    # $ gemini mendel violations
    #########################################
    parser_mendel = subparsers.add_parser('mendel_errors',
            help='Identify candidate violations of Mendelian inheritance')
    parser_mendel.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_mendel.add_argument('--columns',
            dest='columns',
            metavar='STRING',
            help='A list of columns that you would like returned. Def. = "*"',
            )
    parser_mendel.add_argument('--filter',
            dest='filter',
            metavar='STRING',
            help='Restrictions to apply to variants (SQL syntax)')
    parser_mendel.add_argument('-d',
            dest='min_sample_depth',
            type=int,
            help="The minimum aligned\
                  sequence depth (genotype DP) req'd for\
                  each sample (def. = 0)",
            default=0)
    def mendel_fn(parser, args):
        import tool_mendel_errors
        tool_mendel_errors.run(parser, args)
    parser_mendel.set_defaults(func=mendel_fn)


    #########################################
    # $ gemini browser
    #########################################
    parser_browser = subparsers.add_parser('browser',
            help='Browser interface to gemini')
    parser_browser.add_argument('db', metavar='db',
            help='The name of the database to be queried.')
    def browser_fn(parser, args):
        import gemini_browser
        gemini_browser.browser_main(parser, args)
    parser_browser.set_defaults(func=browser_fn)


    #########################################
    # $ gemini set_somatic
    #########################################
    parser_set_somatic = subparsers.add_parser("set_somatic",
                          help="Tag somatic mutations (is_somatic) by comparint tumor/normal pairs.")
    parser_set_somatic.add_argument('db', metavar='db',
            help='The name of the database to be updated.')

    parser_set_somatic.add_argument('--min-depth',
            dest='min_depth',
            type=float,
            default=30,
            help='The min combined depth for tumor + normal (def: %(default)s).')

    parser_set_somatic.add_argument('--min-qual',
            dest='min_qual',
            type=float,
            default=30,
            help='The min variant quality (VCF QUAL) (def: %(default)s).')

    parser_set_somatic.add_argument('--max-norm-alt-freq',
            dest='max_norm_alt_freq',
            type=float,
            default=0.03,
            help='The max freq. of the alt. allele in the normal sample (def: %(default)s).')

    parser_set_somatic.add_argument('--max-norm-alt-count',
            dest='max_norm_alt_count',
            type=int,
            default=2,
            help='The max count. of the alt. allele in the normal sample (def: %(default)s).')

    parser_set_somatic.add_argument('--min-norm-depth',
            dest='min_norm_depth',
            type=int,
            default=10,
            help='The minimum depth allowed in the normal sample to believe somatic (def: %(default)s).')

    parser_set_somatic.add_argument('--min-tumor-alt-freq',
            dest='min_tumor_alt_freq',
            type=float,
            default=0.05,
            help='The min freq. of the alt. allele in the tumor sample (def: %(default)s).')

    parser_set_somatic.add_argument('--min-tumor-alt-count',
            dest='min_tumor_alt_count',
            type=int,
            default=2,
            help='The min count. of the alt. allele in the tumor sample (def: %(default)s).')

    parser_set_somatic.add_argument('--min-tumor-depth',
            dest='min_tumor_depth',
            type=int,
            default=10,
            help='The minimum depth allowed in the tumor sample to believe somatic (def: %(default)s).')

    parser_set_somatic.add_argument('--chrom',
            dest='chrom',
            metavar='STRING',
            help='A specific chromosome on which to tag somatic mutations. (def: %(default)s).',
            default=None,
            )

    parser_set_somatic.add_argument('--dry-run',
            dest='dry_run',
            action='store_true',
            help='Don\'t set the is_somatic flag, just report what _would_ be set. For testing parameters.',
            default=False)

    def set_somatic_fn(parser, args):
        import gemini_set_somatic
        gemini_set_somatic.set_somatic(parser, args)
    parser_set_somatic.set_defaults(func=set_somatic_fn)

    #########################################
    # $ gemini actionable_mutations
    #########################################
    parser_actionable_mut = subparsers.add_parser("actionable_mutations",
                          help="Retrieve genes with actionable somatic mutations via COSMIC and DGIdb.")
    parser_actionable_mut.add_argument('db', metavar='db',
            help='The name of the database to be queried.')
    def get_actionable_mut_fn(parser, args):
        import gemini_actionable_mutations
        gemini_actionable_mutations.get_actionable_mutations(parser, args)
    parser_actionable_mut.set_defaults(func=get_actionable_mut_fn)


    #########################################
    # $ gemini update
    #########################################
    parser_update = subparsers.add_parser("update", help="Update gemini software and data files.")
    parser_update.add_argument("--devel", help="Get the latest development version instead of the release",
                               action="store_true", default=False)
    parser_update.add_argument("--dataonly", help="Only update data, not the underlying libraries.",
                               action="store_true", default=False)
    parser_update.add_argument("--extra", help="Add additional non-standard genome annotations to include",
                               action="append", default=[], choices=["gerp_bp","cadd_score"])
    def update_fn(parser, args):
        import gemini_update
        gemini_update.release(parser, args)
    parser_update.set_defaults(func=update_fn)


    #########################################
    # $ gemini roh
    #########################################
    parser_hom_run = subparsers.add_parser('roh',
            help='Identify runs of homozygosity')
    parser_hom_run.add_argument('db',
            metavar='db',
            help='The name of the database to be queried.')
    parser_hom_run.add_argument('--min-snps',
            dest='min_snps',
            metavar="INTEGER",
            type=int,
            default=25,
            help='Minimum number of homozygous snps expected in a run (def. 25)')
    parser_hom_run.add_argument('--min-total-depth',
            dest='min_total_depth',
            metavar="INTEGER",
            type=int,
            default=20,
            help="""The minimum overall sequencing depth required"""
                 """for a SNP to be considered (def = 20).""")
    parser_hom_run.add_argument('--min-gt-depth',
            dest='min_genotype_depth',
            metavar="INTEGER",
            type=int,
            default=0,
            help="""The minimum required sequencing depth underlying a given sample's genotype"""
                 """for a SNP to be considered (def = 0).""")
    parser_hom_run.add_argument('--min-size',
            metavar="INTEGER",
            dest='min_size',
            type=int,
            default=100000,
            help='Minimum run size in base pairs (def. 100000)')
    parser_hom_run.add_argument('--max-hets',
            metavar="INTEGER",
            dest='max_hets',
            type=int,
            default=1,
            help='Maximum number of allowed hets in the run (def. 1)')
    parser_hom_run.add_argument('--max-unknowns',
            metavar="INTEGER",
            type=int,
            dest='max_unknowns',
            default=3,
            help='Maximum number of allowed unknowns in the run (def. 3)')
    parser_hom_run.add_argument('-s',
            dest='samples',
            default=None,
            help='Comma separated list of samples to screen for ROHs. e.g S120,S450')
    def homozygosity_runs_fn(parser, args):
        from tool_homozygosity_runs import run
        run(parser, args)
    parser_hom_run.set_defaults(func=homozygosity_runs_fn)
    #######################################################
    # parse the args and call the selected function
    #######################################################
    args = parser.parse_args()

    # make sure database is found if provided
    if len(sys.argv) > 2 and sys.argv[1] not in \
       ["load", "merge_chunks", "load_chunk"]:
        if hasattr(args, "db") and args.db is not None and not os.path.exists(args.db):
            sys.stderr.write("Requested GEMINI database (%s) not found. "
                             "Please confirm the provided filename.\n"
                             % args.db)

    try:
        args.func(parser, args)
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

class IPythonAction(argparse.Action):
    def __call__(self, parser, args, values, option = None):
        args.scheduler = values
        if xor(args.scheduler, args.queue):
            parser.error("If you are using the IPython parallel loading, you "
                         "must specify both a scheduler with --scheduler and a "
                         "queue to use with --queue.")

def xor(arg1, arg2):
    return bool(arg1) ^ bool(arg2)


if __name__ == "__main__":
    main()
