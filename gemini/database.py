#!/usr/bin/env python

import os
import contextlib
import sys

import sqlalchemy as sql
from sqlalchemy.orm import mapper, create_session
import sqlalchemy

from ped import get_ped_fields


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
    cursor.execute('''create index chrom_varid_idx on variants(chrom,variant_id)''')
    cursor.execute('CREATE index max_aaf_all_idx on variants(max_aaf_all)')

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


def get_path(path):
    if "://" not in path:
        if not (os.path.exists(path) or os.access(os.path.dirname(os.path.abspath(path)), os.W_OK)):
            raise Exception("unable to find path: %s" % path)
        path = "sqlite:///%s" % path

    return path

def create_tables(path, effect_fields=None):
    """
    Create our master DB tables
    """
    if effect_fields:
        effect_string = "".join(e + " TEXT,\n" for e in effect_fields)
    else:
        effect_string = ""

    db = dict(variants="""
    chrom varchar(15),
    start integer,
    end integer,
    vcf_id text,
    variant_id integer,
    anno_id integer,
    ref text,
    alt text,
    qual float,
    filter text,
    type varchar(15),
    sub_type text,
    gts blob,
    gt_types blob,
    gt_phases blob,
    gt_depths blob,
    gt_ref_depths blob,
    gt_alt_depths blob,
    gt_quals blob,
    gt_copy_numbers blob,
    gt_phred_ll_homref blob,
    gt_phred_ll_het blob,
    gt_phred_ll_homalt blob,
    call_rate float,
    max_aaf_all float,
    in_dbsnp bool,
    rs_ids text default NULL,
    sv_cipos_start_left integer,
    sv_cipos_end_left integer,
    sv_cipos_start_right integer,
    sv_cipos_end_right integer,
    sv_length integer,
    sv_is_precise bool,
    sv_tool text,
    sv_evidence_type text,
    sv_event_id text,
    sv_mate_id text,
    sv_strand text,
    in_omim bool,
    clinvar_sig text default NULL,
    clinvar_disease_name text default NULL,
    clinvar_dbsource text default NULL,
    clinvar_dbsource_id text default NULL,
    clinvar_origin text default NULL,
    clinvar_dsdb text default NULL,
    clinvar_dsdbid text default NULL,
    clinvar_disease_acc text default NULL,
    clinvar_in_locus_spec_db bool,
    clinvar_on_diag_assay bool,
    clinvar_causal_allele text,
    clinvar_gene_phenotype text,
    geno2mp_hpo_ct integer,
    pfam_domain text,
    cyto_band text default NULL,
    rmsk text default NULL,
    in_cpg_island bool,
    in_segdup bool,
    is_conserved bool,
    gerp_bp_score float,
    gerp_element_pval float,
    num_hom_ref integer,
    num_het integer,
    num_hom_alt integer,
    num_unknown integer,
    aaf real,
    hwe decimal(2,7),
    inbreeding_coeff decimal(2,7),
    pi decimal(2,7),
    recomb_rate decimal(2,7),
    gene varchar(60),
    transcript varchar(60),
    is_exonic bool,
    is_coding bool,
    is_splicing bool,
    is_lof bool,
    exon text,
    codon_change text,
    aa_change text,
    aa_length text,
    biotype text,
    impact varchar(60) default NULL,
    impact_so text default NULL,
    impact_severity varchar(15),
    polyphen_pred text,
    polyphen_score float,
    sift_pred text,
    sift_score float,
    anc_allele text,
    rms_bq float,
    cigar text,
    depth integer default NULL,
    strand_bias float default NULL,
    rms_map_qual float default NULL,
    in_hom_run integer default NULL,
    num_mapq_zero integer default NULL,
    num_alleles integer default NULL,
    num_reads_w_dels float default NULL,
    haplotype_score float default NULL,
    qual_depth float default NULL,
    allele_count integer default NULL,
    allele_bal float default NULL,
    in_hm2 bool,
    in_hm3 bool,
    is_somatic bool,
    somatic_score float,
    in_esp bool,
    aaf_esp_ea decimal(2,7),
    aaf_esp_aa decimal(2,7),
    aaf_esp_all decimal(2,7),
    exome_chip bool,
    in_1kg bool,
    aaf_1kg_amr decimal(2,7),
    aaf_1kg_eas decimal(2,7),
    aaf_1kg_sas decimal(2,7),
    aaf_1kg_afr decimal(2,7),
    aaf_1kg_eur decimal(2,7),
    aaf_1kg_all decimal(2,7),
    grc text default NULL,
    gms_illumina float,
    gms_solid float,
    gms_iontorrent float,
    in_cse bool,
    encode_tfbs text,
    encode_dnaseI_cell_count integer,
    encode_dnaseI_cell_list text,
    encode_consensus_gm12878 text,
    encode_consensus_h1hesc text,
    encode_consensus_helas3 text,
    encode_consensus_hepg2 text,
    encode_consensus_huvec text,
    encode_consensus_k562 text,
    vista_enhancers text,
    cosmic_ids text,
    info blob,
    cadd_raw float,
    cadd_scaled float,
    fitcons float,
    in_exac bool,
    aaf_exac_all decimal(2,7),
    aaf_adj_exac_all decimal(2,7),
    aaf_adj_exac_afr decimal(2,7),
    aaf_adj_exac_amr decimal(2,7),
    aaf_adj_exac_eas decimal(2,7),
    aaf_adj_exac_fin decimal(2,7),
    aaf_adj_exac_nfe decimal(2,7),
    aaf_adj_exac_oth decimal(2,7),
    aaf_adj_exac_sas decimal(2,7),
    exac_num_het int,
    exac_num_hom_alt int,
    exac_num_chroms int,
    %s""" % effect_string.rstrip(","),

    variant_impacts="""
    variant_id integer,
    anno_id integer,
    gene varchar(60),
    transcript varchar(60),
    is_exonic bool,
    is_coding bool,
    is_lof bool,
    exon text,
    codon_change text,
    aa_change text,
    aa_length text,
    biotype text,
    impact varchar(60),
    impact_so text,
    impact_severity varchar(15),
    polyphen_pred text,
    polyphen_score float,
    sift_pred text,
    sift_score float""",
    sample_genotypes="""
    sample_id integer,
    gt_types BLOB""",

    sample_genotype_counts="""
    sample_id integer,
    num_hom_ref integer,
    num_het integer,
    num_hom_alt integer,
    num_unknown integer""",

    resources="""
    name text,
    resource text""",

    version="""version text""",

    gene_detailed="""
    uid integer,
    chrom varchar(60),
    gene varchar(60),
    is_hgnc bool,
    ensembl_gene_id text,
    transcript varchar(60),
    biotype text,
    transcript_status text,
    ccds_id varchar(60),
    hgnc_id text,
    entrez_id text,
    cds_length text,
    protein_length text,
    transcript_start text,
    transcript_end text,
    strand text,
    synonym text,
    rvis_pct float,
    mam_phenotype_id text""",

    gene_summary="""
    uid integer,
    chrom varchar(60),
    gene varchar(60),
    is_hgnc bool,
    ensembl_gene_id text,
    hgnc_id text,
    transcript_min_start integer,
    transcript_max_end integer,
    strand text,
    synonym text,
    rvis_pct float,
    mam_phenotype_id text,
    in_cosmic_census bool,
    """,

    vcf_header="""vcf_header text""")

    # in the future this will be replaced by reading from the conf file.

    lookup = {'real': sql.Float(),
              'float': sql.Float(),
              'text': sql.Text(),
              'bool': sql.Boolean(),
              'blob': sql.LargeBinary(),
              'decimal(2,7)': sql.Float(), #sql.DECIMAL(precision=7, scale=2, asdecimal=False),
              'integer': sql.Integer(),
              'varchar(15)': sql.String(20),
              'varchar(60)': sql.String(60),
              'int': sql.Integer(),
              }


    for table in db:

        db[table] = db[table].strip().strip(",").split(",\n")
        db[table] = [x.strip().split() for x in db[table]]
        cols = [sql.Column(c[0], lookup[c[1].lower()]) for c in db[table]]

        if table != "variant_impacts":
            for c in cols:
                if c.name in ("variant_id", "sample_id", "uid"):
                    c.primary_key = True
                    if c.name == "variant_id" and table == "variants":
                        c.autoincrement = False

        db[table] = cols

    e = sql.create_engine(get_path(path), isolation_level=None)
    e.connect().connection.connection.text_factory = str
    metadata = sql.MetaData(bind=e)

    session = create_session(bind=e, autocommit=False, autoflush=False)
    mapped = {}
    tables = ['variants'] + [x for x in sorted(db) if x != 'variants']
    otables = [sql.Table(tbl, metadata, *db[tbl]) for tbl in tables]
    metadata.drop_all(tables=otables)
    for t in otables:
        mapped[t.name] = t
    session.commit()

    metadata.create_all()
    return session, metadata

def create_sample_table(cursor, metadata, args):
    NUM_BUILT_IN = 6
    cols = [sql.Column("sample_id", sql.Integer, primary_key=True)]
    fields = get_ped_fields(args.ped_file)
    optional_fields = ["family_id", "name", "paternal_id", "maternal_id",
                       "sex", "phenotype"]
    optional_fields += fields[NUM_BUILT_IN:]
    for field in optional_fields:
        if field == "name":
            cols.append(sql.Column(field, sql.String(50)))
        else:
            cols.append(sql.Column(field, sql.TEXT))

    t = sql.Table("samples", metadata, *cols)
    t.drop(checkfirst=True)
    metadata.create_all(tables=[t])

def insert_variation(session, metadata, buffer):
    """
    Populate the variants table with each variant in the buffer.
    """
    if len(buffer) == 0: return
    tbl = metadata.tables['variants']

    cols = _get_cols(tbl)

    left = set(cols) - set(buffer[0].keys())
    assert len(left) == 0, left
    left = set(buffer[0].keys()) - set(cols)
    assert len(left) == 0, left

    try:
        session.execute(tbl.insert(), buffer)
        session.commit()
    except:
        sys.stderr.write("insert error trying 1 at a time:\n")
        session.rollback()
        stmt = tbl.insert()
        with session.bind.begin() as trans:
            for b in buffer:
                trans.execute(stmt, b)
        raise


def insert_variation_impacts(session, metadata, buffer):
    """
    Populate the variant_impacts table with each variant in the buffer.
    """
    if len(buffer) == 0: return
    tbl = metadata.tables['variant_impacts']
    session.execute(tbl.insert(), buffer)
    session.commit()


def insert_sample(session, metadata, sample_list):
    """
    Populate the samples with sample ids, names, and
    other indicative information.
    """
    tbl = metadata.tables['samples']
    cols = [c.key for c in tbl.c]
    session.execute(tbl.insert(), dict(zip(cols, sample_list)))


def _get_cols(tbl):
    return [c.key for c in tbl.c]

def gen_gene_vals(cols, table_contents):
    for row in table_contents:
        d = dict(zip(cols, row))
        d['is_hgnc'] = bool(d['is_hgnc'] != '0')
        if d['rvis_pct'] in (None, 'None'):
            d['rvis_pct'] = None
        else:
            d['rvis_pct'] = float(d['rvis_pct'])
        if "in_cosmic_census" in d:
            d['in_cosmic_census'] = bool(d['in_cosmic_census'])

        yield d

def insert_gene_detailed(session, metadata, table_contents):
    t = metadata.tables['gene_detailed']
    cols = _get_cols(t)

    n = 1000
    for chunk in (table_contents[i:i+n] for i in xrange(0, len(table_contents), n)):
        session.execute(t.insert(), list(gen_gene_vals(cols, chunk)))

def insert_gene_summary(session, metadata, contents):
    t = metadata.tables['gene_summary']
    cols = _get_cols(t)

    session.execute(t.insert(), list(gen_gene_vals(cols, contents)))
    session.commit()

def insert_resources(session, metadata, resources):
    """Populate table of annotation resources used in this database.
    """
    t = metadata.tables['resources']
    cols = _get_cols(t)

    session.execute(t.insert(), [dict(zip(cols, r)) for r in resources])
    session.commit()

def insert_vcf_header(session, metadata, vcf_header):
    """Populate a table storing the original VCF header.
    """
    t = metadata.tables['vcf_header']
    session.execute(t.insert(), dict(vcf_header=vcf_header.rstrip("\r\n")))
    session.commit()


def insert_version(session, metadata, version):
    """Populate table of documenting which
    gemini version was used for this database.
    """
    t = metadata.tables['version']
    session.execute(t.insert(), dict(version=version))
    session.commit()


def close_and_commit(session):
    """
    Commit changes to the DB and close out DB cursor.
    """
    session.commit()
    session.close()

def empty_tables(cursor):
    cursor.execute('''delete * from variation''')
    cursor.execute('''delete * from samples''')

def update_gene_summary_w_cancer_census(session, metadata, genes):
    tbl = metadata.tables['gene_summary']
    #http://docs.sqlalchemy.org/en/rel_1_0/core/tutorial.html?highlight=select#inserts-updates-and-deletes
    stmt = tbl.update().where(sql.and_(tbl.c.gene == sql.bindparam("gene_"),
                                   tbl.c.chrom == sql.bindparam("chrom_"))
                              ).values(in_cosmic_census=sql.bindparam("in_cosmic_census_"))
    # need the trailing underscore for bindparam
    cols = ("in_cosmic_census_", "gene_", "chrom_")
    session.execute(stmt, [dict(zip(cols, g)) for g in genes])


def get_session_metadata(path):
    """return an engine"""
    engine = sqlalchemy.create_engine(get_path(path), isolation_level=None)
    engine.connect().connection.connection.text_factory = str
    engine.raw_connection().connection.text_factory = str
    metadata = sql.MetaData(bind=engine)
    metadata.reflect(bind=engine)
    session = create_session(bind=engine, autocommit=False, autoflush=False)
    session.bind.connect().connection.connection.text_factory = str
    session.bind.raw_connection().connection.text_factory = str
    return session, metadata


@contextlib.contextmanager
def database_transaction(db):
    session, metadata = get_session_metadata(db)

    conn = session.connection()
    if session.bind.name == "sqlite":
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode=MEMORY')
    yield conn
    session.commit()
    conn.close()
