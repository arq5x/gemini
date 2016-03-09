#!/usr/bin/env python
import os
import shutil
import sqlite3
import sys
import uuid

import database as gemini_db
import gemini_load_chunk


def append_variant_info(main_curr, chunk_db):
    """
    Append the variant and variant_info data from a chunk_db
    to the main database.
    """

    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    main_curr.execute("BEGIN TRANSACTION")
    cmd = "INSERT INTO variants SELECT * FROM toMerge.variants"
    main_curr.execute(cmd)

    cmd = \
        "INSERT INTO variant_impacts SELECT * FROM toMerge.variant_impacts"
    main_curr.execute(cmd)
    main_curr.execute("END TRANSACTION")

    cmd = "detach toMerge"
    main_curr.execute(cmd)


def append_sample_genotype_counts(main_curr, chunk_db):
    """
    Append the sample_genotype_counts from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "INSERT INTO sample_genotype_counts \
           SELECT * FROM toMerge.sample_genotype_counts"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)


def append_sample_info(main_curr, chunk_db):
    """
    Append the sample info from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "create table samples as select * from toMerge.samples where 1=0"
    main_curr.execute(cmd)

    cmd = "INSERT INTO samples SELECT * FROM toMerge.samples"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)


def append_resource_info(main_curr, chunk_db):
    """
    Append the resource info from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "INSERT INTO resources SELECT * FROM toMerge.resources"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)


def append_version_info(main_curr, chunk_db):
    """
    Append the version info from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "INSERT INTO version SELECT * FROM toMerge.version"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)

def append_vcf_header(main_curr, chunk_db):
    """
    Append the vcf_header from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "INSERT INTO vcf_header SELECT * FROM toMerge.vcf_header"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)

def append_gene_summary(main_curr, chunk_db):
    """
    Append the gene_summary from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "INSERT INTO gene_summary SELECT * FROM toMerge.gene_summary"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)

def append_gene_detailed(main_curr, chunk_db):
    """
    Append the gene_detailed from a chunk_db
    to the main database.
    """
    cmd = "attach ? as toMerge"
    main_curr.execute(cmd, (chunk_db, ))

    cmd = "INSERT INTO gene_detailed SELECT * FROM toMerge.gene_detailed"
    main_curr.execute(cmd)

    cmd = "detach toMerge"
    main_curr.execute(cmd)

def update_sample_genotype_counts(main_curr, chunk_db):
    """
    Update the main sample_genotype_counts table with the
    counts observed in one of the chunked databases (chunk_db)
    """
    curr_db_conn = sqlite3.connect(chunk_db)
    curr_db_conn.isolation_level = None
    curr_db_conn.row_factory = sqlite3.Row
    curr_db_curr = curr_db_conn.cursor()

    cmd = "SELECT sample_id, num_hom_ref, \
                  num_het, num_hom_alt, \
                  num_unknown FROM sample_genotype_counts"
    curr_db_curr.execute(cmd)
    for row in curr_db_curr:
        main_curr.execute("""UPDATE sample_genotype_counts
                          SET num_hom_ref = num_hom_ref + ?,
                              num_het = num_het + ?,
                              num_hom_alt = num_hom_alt + ?,
                              num_unknown = num_unknown + ?
                          WHERE sample_id= ? """,
                          (row['num_hom_ref'],
                           row['num_het'],
                           row['num_hom_alt'],
                           row['num_unknown'],
                           row['sample_id']))
    curr_db_curr.close()


def merge_db_chunks(args):

    # open up a new database
    if os.path.exists(args.db):
        os.remove(args.db)

    gemini_db.create_tables(args.db, gemini_load_chunk.get_extra_effects_fields(args) if args.vcf else [])

    main_conn = sqlite3.connect(args.db)
    main_conn.isolation_level = None
    main_curr = main_conn.cursor()
    main_curr.execute('PRAGMA synchronous = OFF')
    main_curr.execute('PRAGMA journal_mode=MEMORY')

    databases = []
    for database in args.chunkdbs:
        databases.append(database)

    for idx, database in enumerate(databases):

        db = database[0]

        append_variant_info(main_curr, db)

        # we only need to add these tables from one of the chunks.
        if idx == 0:
            append_sample_genotype_counts(main_curr, db)
            append_sample_info(main_curr, db)
            append_resource_info(main_curr, db)
            append_version_info(main_curr, db)
            append_vcf_header(main_curr, db)
            append_gene_summary(main_curr, db)
            append_gene_detailed(main_curr, db)
        else:
            update_sample_genotype_counts(main_curr, db)

    if args.index:
        gemini_db.create_indices(main_curr)

    main_conn.commit()
    main_curr.close()


def merge_chunks(parser, args):
    errors = []
    for try_count in range(2):
        try:
            if try_count > 0:
                tmp_dbs = [os.path.join(args.tempdir, "%s.db" % uuid.uuid4())
                           for _ in args.chunkdbs]
                for chunk_db, tmp_db in zip(args.chunkdbs, tmp_dbs):
                    shutil.copyfile(chunk_db[0], tmp_db)
                    chunk_db[0] = tmp_db

                output_db = args.db
                args.db = os.path.join(args.tempdir, "%s.db" % uuid.uuid4())

            merge_db_chunks(args)

            if try_count > 0:
                shutil.move(args.db, output_db)
                for tmp_db in tmp_dbs:
                    os.remove(tmp_db)
            break
        except sqlite3.OperationalError, e:
            errors.append(str(e))
            sys.stderr.write("sqlite3.OperationalError: %s\n" % e)
    else:
        raise Exception("Attempted workaround for SQLite locking issue on NFS "
                        "drives has failed. One possible reason is that the temp directory "
                        "%s is also on an NFS drive. Error messages from SQLite: %s"
                        % (args.tempdir, " ".join(errors)))
