#!/usr/bin/env python

import sqlite3
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
    cursor.execute('''create index var_in_call_rate_idx on \
                      variants(call_rate)''')
    cursor.execute('''create index var_exonic_idx on variants(is_exonic)''')
    cursor.execute('''create index var_coding_idx on variants(is_coding)''')
    cursor.execute('''create index var_lof_idx on variants(is_lof)''')
    cursor.execute('''create index var_depth_idx on variants(depth)''')
    cursor.execute('''create index var_gene_idx on variants(gene)''')


def index_variation_impacts(cursor):
    cursor.execute('''create index varimp_exonic_idx on \
                      variant_impacts(is_exonic)''')
    cursor.execute('''create index varimp_coding_idx on \
                      variant_impacts(is_coding)''')
    cursor.execute(
        '''create index varimp_lof_idx on variant_impacts(is_lof)''')
    cursor.execute('''create index varimp_impact_idx on \
                      variant_impacts(impact)''')


def index_samples(cursor):
    cursor.execute('''create unique index sample_name_idx on samples(name)''')


def create_indices(cursor):
    """
    Index our master DB tables for speed
    """
    index_variation(cursor)
    index_samples(cursor)


def create_tables(cursor):
    """
    Create our master DB tables
    """
    cursor.execute('''create table if not exists variants  (    \
                    chrom text,                                 \
                    start integer,                              \
                    end integer,                                \
                    vcf_id text,                                \
                    variant_id integer,                         \
                    anno_id integer,                            \
                    ref text,                                   \
                    alt text,                                   \
                    qual float,                                 \
                    filter text,                                \
                    type text,                                  \
                    sub_type text,                              \
                    gts blob,                                   \
                    gt_types blob,                              \
                    gt_phases blob,                             \
                    gt_depths blob,                             \
                    gt_ref_depths blob,                         \
                    gt_alt_depths blob,                         \
                    gt_quals blob,                              \
                    call_rate float,                            \
                    in_dbsnp bool,                              \
                    rs_ids text default NULL,                   \
                    in_omim bool,                               \
                    clinvar_sig text default NULL,              \
                    clinvar_disease_name text default NULL,     \
                    clinvar_dbsource text default NULL,         \
                    clinvar_dbsource_id text default NULL,      \
                    clinvar_origin text default NULL,           \
                    clinvar_dsdb text default NULL,             \
                    clinvar_dsdbid text default NULL,           \
                    clinvar_disease_acc text default NULL,      \
                    clinvar_in_locus_spec_db bool,              \
                    clinvar_on_diag_assay bool,                 \
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
                    hwe decimal(2,7),                           \
                    inbreeding_coeff decimal(2,7),              \
                    pi decimal(2,7),                            \
                    recomb_rate decimal(2,7),                   \
                    gene text,                                  \
                    transcript text,                            \
                    is_exonic bool,                             \
                    is_coding bool,                             \
                    is_lof bool,                                \
                    exon text,                                  \
                    codon_change text,                          \
                    aa_change text,                             \
                    aa_length text,                             \
                    biotype text,                               \
                    impact text default NULL,                   \
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
                    in_hm2 bool,                                \
                    in_hm3 bool,                                \
                    is_somatic,                                 \
                    in_esp bool,                                \
                    aaf_esp_ea decimal(2,7),                    \
                    aaf_esp_aa decimal(2,7),                    \
                    aaf_esp_all decimal(2,7),                   \
                    exome_chip bool,                            \
                    in_1kg bool,                                \
                    aaf_1kg_amr decimal(2,7),                   \
                    aaf_1kg_asn decimal(2,7),                   \
                    aaf_1kg_afr decimal(2,7),                   \
                    aaf_1kg_eur decimal(2,7),                   \
                    aaf_1kg_all decimal(2,7),                   \
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
                    PRIMARY KEY(variant_id ASC))''')

    cursor.execute('''create table if not exists variant_impacts  (   \
                    variant_id integer,                               \
                    anno_id integer,                                  \
                    gene text,                                        \
                    transcript text,                                  \
                    is_exonic bool,                                   \
                    is_coding bool,                                   \
                    is_lof bool,                                      \
                    exon text,                                        \
                    codon_change text,                                \
                    aa_change text,                                   \
                    aa_length text,                                   \
                    biotype text,                                     \
                    impact text,                                      \
                    impact_severity text,                             \
                    polyphen_pred text,                               \
                    polyphen_score float,                             \
                    sift_pred text,                                   \
                    sift_score float,                                 \
                    PRIMARY KEY(variant_id ASC, anno_id ASC))''')

    cursor.execute('''create table if not exists sample_genotypes (  \
                    sample_id integer,                               \
                    gt_types BLOB,                                   \
                    PRIMARY KEY(sample_id ASC))''')

    cursor.execute('''create table if not exists sample_genotype_counts ( \
                     sample_id integer,                                   \
                     num_hom_ref integer,                                 \
                     num_het integer,                                     \
                     num_hom_alt integer,                                 \
                     num_unknown integer,                                 \
                     PRIMARY KEY(sample_id ASC))''')

    cursor.execute('''create table if not exists resources ( \
                     name text,                              \
                     resource text)''')

    cursor.execute('''create table if not exists version (version text)''')

def create_sample_table(cursor, args):
    NUM_BUILT_IN = 6
    fields = get_ped_fields(args.ped_file)
    required = "sample_id integer"
    optional_fields = ["family_id", "name", "paternal_id", "maternal_id",
                       "sex", "phenotype"]
    optional_fields += fields[NUM_BUILT_IN:] + ["PRIMARY KEY(sample_id ASC)"]
    optional = " text default NULL,".join(optional_fields)
    structure = '''{0}, {1}'''.format(required, optional)
    creation = "create table if not exists samples ({0})".format(structure)
    cursor.execute(creation)

def _insert_variation_one_per_transaction(cursor, buffer):
    for variant in buffer:
        try:
            cursor.execute("BEGIN TRANSACTION")
            cursor.execute('insert into variants values     (?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?,?,?, \
                                                             ?,?,?,?,?,?,?,?)', variant)
            cursor.execute("END TRANSACTION")
        # skip repeated keys until we get to the failed variant
        except sqlite3.IntegrityError, e:
            cursor.execute("END TRANSACTION")
            continue
        except sqlite3.ProgrammingError, e:
            print variant
            print "Error %s:" % (e.args[0])
            sys.exit(1)

def insert_variation(cursor, buffer):
    """
    Populate the variants table with each variant in the buffer.
    """
    try:
        cursor.execute("BEGIN TRANSACTION")
        cursor.executemany('insert into variants values (?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?,?,?, \
                                                         ?,?,?,?,?,?,?,?)', buffer)

        cursor.execute("END TRANSACTION")
    except sqlite3.ProgrammingError:
        cursor.execute("END TRANSACTION")
        _insert_variation_one_per_transaction(cursor, buffer)


def insert_variation_impacts(cursor, buffer):
    """
    Populate the variant_impacts table with each variant in the buffer.
    """
    cursor.execute("BEGIN TRANSACTION")
    cursor.executemany('insert into variant_impacts values (?,?,?,?,?,?,?,?, \
                                                            ?,?,?,?,?,?,?,?, \
                                                            ?,?)',
                       buffer)
    cursor.execute("END")


def insert_sample(cursor, sample_list):
    """
    Populate the samples with sample ids, names, and
    other indicative information.
    """
    placeholders = ",".join(list(repeat("?", len(sample_list))))
    cursor.execute("BEGIN TRANSACTION")
    cursor.execute("insert into samples values "
                   "({0})".format(placeholders), sample_list)
    cursor.execute("END")


def insert_resources(cursor, resources):
    """Populate table of annotation resources used in this database.
    """
    cursor.execute("BEGIN TRANSACTION")
    cursor.executemany('''insert into resources values (?,?)''', resources)
    cursor.execute("END")

def insert_version(cursor, version):
    """Populate table of documenting which
    gemini version was used for this database.
    """
    cursor.execute("BEGIN TRANSACTION")
    cursor.execute('''insert into version values (?)''', (version,))
    cursor.execute("END")


def close_and_commit(cursor, connection):
    """
    Commit changes to the DB and close out DB cursor.
    """
    connection.commit()
    cursor.close


def empty_tables(cursor):
    cursor.execute('''delete * from variation''')
    cursor.execute('''delete * from samples''')


@contextlib.contextmanager
def database_transaction(db):
    conn = sqlite3.connect(db)
    conn.isolation_level = None
    cursor = conn.cursor()
    cursor.execute('PRAGMA synchronous = OFF')
    cursor.execute('PRAGMA journal_mode=MEMORY')
    yield cursor
    conn.commit
    cursor.close()
