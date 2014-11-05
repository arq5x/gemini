#!/usr/bin/env python

import  psycopg2
import sys
from itertools import repeat
import contextlib

from ped import get_ped_fields, default_ped_fields


def index_variation(cursor):
    cursor.execute('''create index var_chr_start_idx on\
                      variants(chrom, start)''')
    cursor.execute('''create index var_type_idx on variants(type)''')
    cursor.execute('''create index var_gt_counts_idx on \
                      variants(num_hom_ref, num_het, \
                               num_hom_alt, num_unknown)''')
    cursor.execute('''create index var_aaf_idx on variants(aaf)''')
    cursor.execute('''create index var_in_dbsnp_idx on variants(in_dbsnp)''')
    cursor.execute('''create index var_in_call_rate_idx on variants(call_rate)''')
    cursor.execute('''create index var_exonic_idx on variants(is_exonic)''')
    cursor.execute('''create index var_coding_idx on variants(is_coding)''')
    cursor.execute('''create index var_lof_idx on variants(is_lof)''')
    cursor.execute('''create index var_som_idx on variants(is_somatic)''')
    cursor.execute('''create index var_depth_idx on variants(depth)''')
    cursor.execute('''create index var_gene_idx on variants(gene)''')
    cursor.execute('''create index var_trans_idx on variants(transcript)''')
    cursor.execute('''create index var_impact_idx on variants(impact)''')
    cursor.execute('''create index var_impact_severity_idx on variants(impact_severity)''')
    cursor.execute('''create index var_esp_idx on variants(aaf_esp_all)''')
    cursor.execute('''create index var_1kg_idx on variants(aaf_1kg_all)''')
    cursor.execute('''create index var_qual_idx on variants(qual)''')
    cursor.execute('''create index var_homref_idx on variants(num_hom_ref)''')
    cursor.execute('''create index var_homalt_idx on variants(num_hom_alt)''')
    cursor.execute('''create index var_het_idx on variants(num_het)''')
    cursor.execute('''create index var_unk_idx on variants(num_unknown)''')
    cursor.execute('''create index var_omim_idx on variants(in_omim)''')
    cursor.execute('''create index var_cadd_raw_idx on variants(cadd_raw)''')
    cursor.execute('''create index var_cadd_scaled_idx on variants(cadd_scaled)''')
    cursor.execute('''create index var_fitcons_idx on variants(fitcons)''')
    cursor.execute('''create index var_sv_event_idx on variants(sv_event_id)''')
    # genotype array indices
    cursor.execute('''create index vid_gt_idx on variants using GIN (gts)''')
    cursor.execute('''create index vid_gt_types_idx on variants using GIN (gt_types)''')


def index_variation_impacts(cursor):
    cursor.execute('''create index varimp_exonic_idx on \
                      variant_impacts(is_exonic)''')
    cursor.execute('''create index varimp_coding_idx on \
                      variant_impacts(is_coding)''')
    cursor.execute(
        '''create index varimp_lof_idx on variant_impacts(is_lof)''')
    cursor.execute('''create index varimp_impact_idx on \
                      variant_impacts(impact)''')
    cursor.execute('''create index varimp_trans_idx on \
                      variant_impacts(transcript)''')
    cursor.execute('''create index varimp_gene_idx on \
                      variant_impacts(gene)''')


def index_samples(cursor):
    cursor.execute('''create unique index sample_name_idx on samples(name)''')


def index_gene_detailed(cursor):
    cursor.execute('''create index gendet_chrom_gene_idx on \
                       gene_detailed(chrom, gene)''')
    cursor.execute('''create index gendet_rvis_idx on \
                       gene_detailed(rvis_pct)''')
    cursor.execute('''create index gendet_transcript_idx on \
                       gene_detailed(transcript)''')
    cursor.execute('''create index gendet_ccds_idx on \
                       gene_detailed(ccds_id)''')

def index_gene_summary(cursor):
    cursor.execute('''create index gensum_chrom_gene_idx on \
                       gene_summary(chrom, gene)''')
    cursor.execute('''create index gensum_rvis_idx on \
                      gene_summary(rvis_pct)''')

def create_indices(cursor):
    """
    Index our master DB tables for speed
    """
    index_variation(cursor)
    index_variation_impacts(cursor)
    index_samples(cursor)
    index_gene_detailed(cursor)
    index_gene_summary(cursor)


def drop_tables(cursor):
    cursor.execute("DROP TABLE IF EXISTS variants")
    cursor.execute("DROP TABLE IF EXISTS variant_impacts")
    cursor.execute("DROP TABLE IF EXISTS sample_genotypes")
    cursor.execute("DROP TABLE IF EXISTS sample_genotype_counts")
    cursor.execute("DROP TABLE IF EXISTS resources")
    cursor.execute("DROP TABLE IF EXISTS version")
    cursor.execute("DROP TABLE IF EXISTS gene_detailed")
    cursor.execute("DROP TABLE IF EXISTS gene_summary")



def create_tables(cursor):
    """
    Create our master DB tables
    """
    cursor.execute('''create table if not exists variants  (                  \
                    chrom text,                                 \
                    start integer,                              \
                    \"end\" integer,                            \
                    vcf_id text,                                \
                    variant_id integer,                         \
                    anno_id integer,                            \
                    ref text,                                   \
                    alt text,                                   \
                    qual float,                                 \
                    filter text,                                \
                    type text,                                  \
                    sub_type text,                              \
                    gts text[],                                 \
                    gt_types integer[],                         \
                    gt_phases integer[],                        \
                    gt_depths integer[],                        \
                    gt_ref_depths integer[],                    \
                    gt_alt_depths integer[],                    \
                    gt_quals float[],                           \
                    gt_copy_numbers float[],                    \
                    call_rate float,                            \
                    in_dbsnp integer,                           \
                    rs_ids text default NULL,                   \
                    sv_cipos_start_left integer,                \
                    sv_cipos_end_left integer,                  \
                    sv_cipos_start_right integer,               \
                    sv_cipos_end_right integer,                 \
                    sv_length integer,                          \
                    sv_is_precise bool,                         \
                    sv_tool text,                               \
                    sv_evidence_type text,                      \
                    sv_event_id text,                           \
                    sv_mate_id text,                            \
                    sv_strand text,                             \
                    in_omim integer,                            \
                    clinvar_sig text default NULL,              \
                    clinvar_disease_name text default NULL,     \
                    clinvar_dbsource text default NULL,         \
                    clinvar_dbsource_id text default NULL,      \
                    clinvar_origin text default NULL,           \
                    clinvar_dsdb text default NULL,             \
                    clinvar_dsdbid text default NULL,           \
                    clinvar_disease_acc text default NULL,      \
                    clinvar_in_locus_spec_db integer,           \
                    clinvar_on_diag_assay integer,              \
                    clinvar_causal_allele text,                 \
                    pfam_domain text,                           \
                    cyto_band text default NULL,                \
                    rmsk text default NULL,                     \
                    in_cpg_island bool,                         \
                    in_segdup bool,                             \
                    is_conserved bool,                          \
                    gerp_bp_score float,                        \
                    gerp_element_pval float,                    \
                    num_hom_ref integer,                        \
                    num_het integer,                            \
                    num_hom_alt integer,                        \
                    num_unknown integer,                        \
                    aaf real,                                   \
                    hwe decimal(9,7),                           \
                    inbreeding_coeff numeric,                   \
                    pi numeric,                                 \
                    recomb_rate numeric,                        \
                    gene text,                                  \
                    transcript text,                            \
                    is_exonic integer,                             \
                    is_coding integer,                             \
                    is_lof integer,                                \
                    exon text,                                  \
                    codon_change text,                          \
                    aa_change text,                             \
                    aa_length text,                             \
                    biotype text,                               \
                    impact text default NULL,                   \
                    impact_so text default NULL,                \
                    impact_severity text,                       \
                    polyphen_pred text,                         \
                    polyphen_score float,                       \
                    sift_pred text,                             \
                    sift_score float,                           \
                    anc_allele text,                            \
                    rms_bq float,                               \
                    cigar text,                                 \
                    depth integer default NULL,                 \
                    strand_bias float default NULL,             \
                    rms_map_qual float default NULL,            \
                    in_hom_run integer default NULL,            \
                    num_mapq_zero integer default NULL,         \
                    num_alleles integer default NULL,           \
                    num_reads_w_dels float default NULL,        \
                    haplotype_score float default NULL,         \
                    qual_depth float default NULL,              \
                    allele_count integer default NULL,          \
                    allele_bal float default NULL,              \
                    in_hm2 integer,                                \
                    in_hm3 integer,                                \
                    is_somatic integer,                            \
                    somatic_score float,                        \
                    in_esp bool,                                \
                    aaf_esp_ea numeric,                    \
                    aaf_esp_aa numeric,                    \
                    aaf_esp_all numeric,                   \
                    exome_chip bool,                            \
                    in_1kg bool,                                \
                    aaf_1kg_amr numeric,                   \
                    aaf_1kg_eas numeric,                   \
                    aaf_1kg_sas numeric,                   \
                    aaf_1kg_afr numeric,                   \
                    aaf_1kg_eur numeric,                   \
                    aaf_1kg_all numeric,                   \
                    grc text default NULL,                      \
                    gms_illumina float,                         \
                    gms_solid float,                            \
                    gms_iontorrent float,                       \
                    in_cse bool,                                \
                    encode_tfbs text,                           \
                    encode_dnaseI_cell_count integer,           \
                    encode_dnaseI_cell_list text,               \
                    encode_consensus_gm12878 text,              \
                    encode_consensus_h1hesc text,               \
                    encode_consensus_helas3 text,               \
                    encode_consensus_hepg2 text,                \
                    encode_consensus_huvec text,                \
                    encode_consensus_k562 text,                 \
                    vista_enhancers text,                       \
                    cosmic_ids text,                            \
                    info BYTEA,                                  \
                    cadd_raw float,                             \
                    cadd_scaled float,                          \
                    fitcons float,                              \
                    PRIMARY KEY(variant_id))''')

    cursor.execute('''create table if not exists variant_impacts  (   \
                    variant_id integer,                               \
                    anno_id integer,                                  \
                    gene text,                                        \
                    transcript text,                                  \
                    is_exonic integer,                                   \
                    is_coding integer,                                   \
                    is_lof integer,                                      \
                    exon text,                                        \
                    codon_change text,                                \
                    aa_change text,                                   \
                    aa_length text,                                   \
                    biotype text,                                     \
                    impact text,                                      \
                    impact_so text,                                   \
                    impact_severity text,                             \
                    polyphen_pred text,                               \
                    polyphen_score float,                             \
                    sift_pred text,                                   \
                    sift_score float,                                 \
                    PRIMARY KEY(variant_id, anno_id))''')

    cursor.execute('''create table if not exists sample_genotypes (  \
                    sample_id integer,                               \
                    gt_types integer[],                                   \
                    PRIMARY KEY(sample_id))''')

    cursor.execute('''create table if not exists sample_genotype_counts ( \
                     sample_id integer,                                   \
                     num_hom_ref integer,                                 \
                     num_het integer,                                     \
                     num_hom_alt integer,                                 \
                     num_unknown integer,                                 \
                     PRIMARY KEY(sample_id))''')

    cursor.execute('''create table if not exists resources ( \
                     name text,                              \
                     resource text)''')

    cursor.execute('''create table if not exists version (version text)''')
    
    cursor.execute('''create table if not exists gene_detailed (       \
                   uid integer,                                        \
                   chrom text,                                         \
                   gene text,                                          \
                   is_hgnc integer,                                       \
                   ensembl_gene_id text,                               \
                   transcript text,                                    \
                   biotype text,                                       \
                   transcript_status text,                             \
                   ccds_id text,                                       \
                   hgnc_id text,                                       \
                   entrez_id text,                                     \
                   cds_length text,                                    \
                   protein_length text,                                \
                   transcript_start text,                              \
                   transcript_end text,                                \
                   strand text,                                        \
                   synonym text,                                       \
                   rvis_pct float,                                     \
                   mam_phenotype_id text,                              \
                   PRIMARY KEY(uid))''')
    
    cursor.execute('''create table if not exists gene_summary (     \
                    uid integer,                                    \
                    chrom text,                                     \
                    gene text,                                      \
                    is_hgnc integer,                                   \
                    ensembl_gene_id text,                           \
                    hgnc_id text,                                   \
                    transcript_min_start text,                      \
                    transcript_max_end text,                        \
                    strand text,                                    \
                    synonym text,                                   \
                    rvis_pct float,                                 \
                    mam_phenotype_id text,                          \
                    in_cosmic_census integer,                          \
                    PRIMARY KEY(uid))''')

def create_sample_table(cursor, args):
    NUM_BUILT_IN = 6
    fields = get_ped_fields(args.ped_file)
    required = "sample_id integer"
    optional_fields = ["family_id", "name", "paternal_id", "maternal_id",
                       "sex", "phenotype"]
    optional_fields += fields[NUM_BUILT_IN:] + ["PRIMARY KEY(sample_id)"]
    optional = " text default NULL,".join(optional_fields)
    structure = '''{0}, {1}'''.format(required, optional)
    creation = "create table if not exists samples ({0})".format(structure)
    cursor.execute(creation)

def _insert_variation_one_per_transaction(cursor, buffer):
    for variant in buffer:
        try:
            cursor.execute("BEGIN TRANSACTION")
            cursor.execute('insert into variants values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', variant)
            cursor.execute("END TRANSACTION")
        # skip repeated keys until we get to the failed variant
        except psycopg2.IntegrityError, e:
            cursor.execute("END TRANSACTION")
            continue
        except psycopg2.ProgrammingError, e:
            print variant
            print "Error %s:" % (e.args[0])
            sys.exit(1)

def insert_variation(cursor, buffer):
    """
    Populate the variants table with each variant in the buffer.
    """
    try:
        cursor.execute("BEGIN TRANSACTION")
        cursor.executemany('insert into variants values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', buffer)

        cursor.execute("END TRANSACTION")
    except psycopg2.ProgrammingError:
        cursor.execute("END TRANSACTION")
        _insert_variation_one_per_transaction(cursor, buffer)


def insert_variation_impacts(cursor, buffer):
    """
    Populate the variant_impacts table with each variant in the buffer.
    """
    cursor.execute("BEGIN TRANSACTION")
    cursor.executemany('insert into variant_impacts values (%s,%s,%s,%s,%s,%s,%s,%s, \
                                                            %s,%s,%s,%s,%s,%s,%s,%s, \
                                                            %s,%s,%s)',
                       buffer)
    cursor.execute("END")


def insert_sample(cursor, sample_list):
    """
    Populate the samples with sample ids, names, and
    other indicative information.
    """
    placeholders = ",".join(list(repeat("%s", len(sample_list))))
    cursor.execute("BEGIN TRANSACTION")
    
    # a hack to prevent loading the same data multiple times in PGSQL mode.
    cursor.execute("select count(1) from samples")
    count = cursor.fetchone()[0]
    if count == 0:
        cursor.execute("insert into samples values "
                       "({0})".format(placeholders), sample_list)
    
    cursor.execute("END")

def insert_gene_detailed(cursor, table_contents):
    cursor.execute("BEGIN TRANSACTION")
    cursor.executemany('insert into gene_detailed values (%s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                          %s,%s,%s,%s,%s,%s,%s,%s,%s, \
                                                          %s)',
                        table_contents)
    cursor.execute("END")
    

def insert_gene_summary(cursor, contents):
    cursor.execute("BEGIN TRANSACTION")
    cursor.executemany('insert into gene_summary values (%s,%s,%s,%s,%s,%s,%s,%s, \
                                                         %s,%s,%s,%s,%s)', 
                        contents)
    cursor.execute("END")
    
def insert_resources(cursor, resources):
    """Populate table of annotation resources used in this database.
    """
    cursor.execute("BEGIN TRANSACTION")
    cursor.executemany('''insert into resources values (%s,%s)''', resources)
    cursor.execute("END")

def insert_version(cursor, version):
    """Populate table of documenting which
    gemini version was used for this database.
    """
    cursor.execute("BEGIN TRANSACTION")
    cursor.execute('''insert into version values (%s)''', (version,))
    cursor.execute("END")


def close_and_commit(cursor, connection):
    """
    Commit changes to the DB and close out DB cursor.
    """
    print "committing"
    connection.commit()

    #cursor.execute("""SELECT * FROM variants 
    #                        WHERE gt_types[7] =1 
    #                        AND   gt_types[9] =0 
    #                        AND   gt_types[17] =2""")
    #for row in cursor:
    #    print row

    #print "closing"
    connection.close()


def empty_tables(cursor):
    cursor.execute('''delete * from variation''')
    cursor.execute('''delete * from samples''')


def update_gene_summary_w_cancer_census(cursor, genes):
    update_qry = "UPDATE gene_summary SET in_cosmic_census = %s "
    update_qry += " WHERE gene = %s and chrom = %s"
    cursor.executemany(update_qry, genes)

# @contextlib.contextmanager
# def database_transaction(db):
#     conn = sqlite3.connect(db)
#     conn.isolation_level = None
#     cursor = conn.cursor()
#     cursor.execute('PRAGMA synchronous = OFF')
#     cursor.execute('PRAGMA journal_mode=MEMORY')
#     yield cursor
#     conn.commit
#     cursor.close()
