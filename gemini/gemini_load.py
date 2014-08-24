#!/usr/bin/env python

# native Python imports
import os.path
import sys
import sqlite3

import annotations
import subprocess
from cluster_helper.cluster import cluster_view
import database as gemini_db
from gemini_load_chunk import GeminiLoader
import gemini_annotate
import uuid
import time
import datetime


def load(parser, args):
    if (args.db is None or args.vcf is None):
        parser.print_help()
        exit("ERROR: load needs both a VCF file and a database file\n")

    annos = annotations.get_anno_files( args )
    # force skipping CADD and GERP if the data files have not been installed
    if args.skip_cadd is False:
        if 'cadd_score' not in annos:
            sys.stderr.write("\nCADD scores are not being loaded because the"
            " annotation file could not be found.\n"
            "`Run gemini update --dataonly --extra cadd_score`"
            " to install the annotation file.\n\n")
            args.skip_cadd = True
        else:
            sys.stderr.write("CADD scores are being loaded (to skip use:--skip-cadd).\n")
    if args.skip_gerp_bp is False:
        if 'gerp_bp' not in annos:
            sys.stderr.write("\nGERP per bp is not being loaded because the annotation file"
                        " could not be found.\n    Run `gemini update --dataonly --extra gerp_bp`"
                        " to install the annotation file.\n\n")
            args.skip_gerp_bp = True
        else:
            sys.stderr.write("GERP per bp is being loaded (to skip use:--skip-gerp-bp).\n")
    # collect of the the add'l annotation files
    annotations.load_annos( args )

    if args.scheduler:
        load_ipython(args)
    elif args.cores > 1:
        load_multicore(args)
    else:
        load_singlecore(args)

def load_singlecore(args):
    # create a new gemini loader and populate
    # the gemini db and files from the VCF
    gemini_loader = GeminiLoader(args)
    gemini_loader.store_resources()
    gemini_loader.store_version()
    gemini_loader.populate_from_vcf()


    if not args.skip_gene_tables and not args.test_mode:
        gemini_loader.update_gene_table()
    if not args.test_mode:
        gemini_loader.build_indices_and_disconnect()

    if not args.no_genotypes and not args.no_load_genotypes:
        gemini_loader.store_sample_gt_counts()
    gemini_annotate.add_extras(args.db, [args.db])

def load_multicore(args):
    grabix_file = bgzip(args.vcf)
    chunks = load_chunks_multicore(grabix_file, args)
    merge_chunks_multicore(chunks, args.db)
    gemini_annotate.add_extras(args.db, chunks)

def load_ipython(args):
    grabix_file = bgzip(args.vcf)
    with cluster_view(*get_ipython_args(args)) as view:
        chunks = load_chunks_ipython(grabix_file, args, view)
        merge_chunks_ipython(chunks, args.db, view)
    gemini_annotate.add_extras(args.db, chunks)

def merge_chunks(chunks, db):
    cmd = get_merge_chunks_cmd(chunks, db)
    print "Merging chunks."
    subprocess.check_call(cmd, shell=True)
    cleanup_temp_db_files(chunks)
    return db

def get_merge_chunks_cmd(chunks, db):
    chunk_names = ""
    for chunk in chunks:
        chunk_names += " --chunkdb  " + chunk
    return "gemini merge_chunks {chunk_names} --db {db}".format(**locals())

def merge_chunks_ipython(chunks, db, view):
    if len(chunks) == 1:
        os.rename(chunks[0], db)
        return db
    else:
        sub_merges = get_chunks_to_merge(chunks)
        tmp_dbs = get_temp_dbs(len(sub_merges), os.path.dirname(sub_merges[0][0]))
        view.map(merge_chunks, sub_merges, tmp_dbs)
        merge_chunks_ipython(tmp_dbs, db, view)

def merge_chunks_multicore(chunks, db):
    if len(chunks) <= 1:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print st, "indexing final database."

        os.rename(chunks[0], db)
        main_conn = sqlite3.connect(db)
        main_conn.isolation_level = None
        main_curr = main_conn.cursor()
        main_curr.execute('PRAGMA synchronous = OFF')
        main_curr.execute('PRAGMA journal_mode=MEMORY')
        
        gemini_db.create_indices(main_curr)
        
        main_conn.commit()
        main_curr.close()
        return db
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print st, "merging", len(chunks), "chunks."
        procs = []
        sub_merges = get_chunks_to_merge(chunks)
        tmp_dbs = get_temp_dbs(len(sub_merges), os.path.dirname(sub_merges[0][0]))
        for sub_merge, tmp_db in zip(sub_merges, tmp_dbs):
            cmd = get_merge_chunks_cmd(sub_merge, tmp_db)
            procs.append(subprocess.Popen(cmd, shell=True))
        wait_until_finished(procs)
        cleanup_temp_db_files(chunks)
        merge_chunks_multicore(tmp_dbs, db)

def get_chunks_to_merge(chunks):
    sublist = list_to_sublists(chunks, 2)
    if len(sublist[-1]) > 1:
        return sublist
    else:
        sublist[-2].extend(sublist[-1])
        return sublist[:-1]

def list_to_sublists(l, n):
    """ convert list l to sublists of length n """
    return [l[i:i+n] for i in xrange(0, len(l), n)]

def get_temp_dbs(n, tmp_dir):
    return [os.path.join(tmp_dir, str(uuid.uuid4())) + ".db" for _ in xrange(n)]

def get_chunk_name(chunk):
    return "--chunkdb " + chunk

def load_chunks_multicore(grabix_file, args):
    cores = args.cores

    # specify the PED file if given one
    ped_file = ""
    if args.ped_file is not None:
        ped_file = "-p " + args.ped_file

    # specify the annotation type if given one
    anno_type = ""
    if args.anno_type is not None:
        anno_type = "-t " + args.anno_type

    no_genotypes = ""
    if args.no_genotypes is True:
        no_genotypes = "--no-genotypes"

    no_load_genotypes = ""
    if args.no_load_genotypes is True:
        no_load_genotypes = "--no-load-genotypes"

    skip_gerp_bp = ""
    if args.skip_gerp_bp is True:
        skip_gerp_bp = "--skip-gerp-bp"

    skip_gene_tables = ""
    if args.skip_gene_tables is True:
        skip_gene_tables = "--skip-gene-tables"

    skip_cadd = ""
    if args.skip_cadd is True:
        skip_cadd = "--skip-cadd"
    
    test_mode = ""
    if args.test_mode is True:
        test_mode = "--test-mode"

    passonly = ""
    if args.passonly is True:
        passonly = "--passonly"

    skip_info_string = ""
    if args.skip_info_string is True:
        skip_info_string = "--skip-info-string"

    submit_command = "{cmd}"
    vcf, _ = os.path.splitext(grabix_file)
    chunk_steps = get_chunk_steps(grabix_file, args)
    chunk_vcfs = []
    chunk_dbs = []
    procs = []

    for chunk_num, chunk in chunk_steps:
        start, stop = chunk
        print "Loading chunk " + str(chunk_num) + "."
        gemini_load = gemini_pipe_load_cmd().format(**locals())
        procs.append(subprocess.Popen(submit_command.format(cmd=gemini_load),
                                      shell=True))

        chunk_vcf = vcf + ".chunk" + str(chunk_num)
        chunk_vcfs.append(chunk_vcf)
        chunk_dbs.append(chunk_vcf + ".db")

    wait_until_finished(procs)
    print "Done loading {0} variants in {1} chunks.".format(stop, chunk_num+1)
    return chunk_dbs

def load_chunks_ipython(grabix_file, args, view):
    # specify the PED file if given one
    ped_file = ""
    if args.ped_file is not None:
        ped_file = "-p " + args.ped_file

    # specify the annotation type if given one
    anno_type = ""
    if args.anno_type is not None:
        anno_type = "-t " + args.anno_type

    no_genotypes = ""
    if args.no_genotypes is True:
        no_genotypes = "--no-genotypes"

    no_load_genotypes = ""
    if args.no_load_genotypes is True:
        no_load_genotypes = "--no-load-genotypes"

    skip_gerp_bp = ""
    if args.skip_gerp_bp is True:
        skip_gerp_bp = "--skip-gerp-bp"

    skip_gene_tables = ""
    if args.skip_gene_tables is True:
        skip_gene_tables = "--skip-gene-tables"
    
    skip_cadd = ""
    if args.skip_cadd is True:
        skip_cadd = "--skip-cadd"
    
    test_mode = ""
    if args.test_mode is True:
        test_mode = "--test-mode"

    passonly = ""
    if args.passonly is True:
        passonly = "--passonly"

    skip_info_string = ""
    if args.skip_info_string is True:
        skip_info_string = "--skip-info-string"


    vcf, _ = os.path.splitext(grabix_file)
    chunk_steps = get_chunk_steps(grabix_file, args)
    total_chunks = len(chunk_steps)
    scheduler, queue, cores = get_ipython_args(args)
    load_args = {"ped_file": ped_file,
                 "anno_type": anno_type,
                 "vcf": vcf,
                 "grabix_file": grabix_file,
                 "no_genotypes": no_genotypes,
                 "no_load_genotypes": no_load_genotypes,
                 "skip_gerp_bp": skip_gerp_bp,
                 "skip_gene_tables": skip_gene_tables,
                 "skip_cadd": skip_cadd,
                 "test_mode": test_mode,
                 "passonly": passonly,
                 "skip_info_string": skip_info_string}
    chunk_dbs = view.map(load_chunk, chunk_steps, [load_args] * total_chunks)

    print "Done loading variants in {0} chunks.".format(total_chunks)
    return chunk_dbs

def load_chunk(chunk_step, kwargs):
    chunk_num, chunk = chunk_step
    start, stop = chunk
    args = combine_dicts(locals(), kwargs)
    gemini_load = gemini_pipe_load_cmd().format(**args)
    subprocess.check_call(gemini_load, shell=True)
    chunk_db = args["vcf"] + ".chunk" + str(chunk_num) + ".db"
    return chunk_db

def wait_until_finished(procs):
    [p.wait() for p in procs]

def cleanup_temp_db_files(chunk_dbs):
    for chunk_db in chunk_dbs:
        os.remove(chunk_db)

def gemini_pipe_load_cmd():
    grabix_cmd = "grabix grab {grabix_file} {start} {stop}"
    gemini_load_cmd = ("gemini load_chunk -v - {anno_type} {ped_file}"
                       " {no_genotypes} {no_load_genotypes} {no_genotypes}"
                       " {skip_gerp_bp} {skip_gene_tables} {skip_cadd}"
                       " {passonly} {skip_info_string} {test_mode}"
                       " -o {start} {vcf}.chunk{chunk_num}.db")
    return " | ".join([grabix_cmd, gemini_load_cmd])

def get_chunk_steps(grabix_file, args):
    index_file = grabix_index(grabix_file)
    num_lines = get_num_lines(index_file)
    chunk_size = int(num_lines) / int(args.cores)
    print "Breaking {0} into {1} chunks.".format(grabix_file, args.cores)

    starts = []
    stops = []
    for chunk in range(0, int(args.cores)):
        start = (chunk * chunk_size) + 1
        stop  = start + chunk_size - 1
        # make sure the last chunk covers the remaining lines
        if chunk >= (args.cores - 1) and stop < num_lines:
            stop = num_lines
        starts.append(start)
        stops.append(stop)
    return list(enumerate(zip(starts, stops)))

def get_num_lines(index_file):
    with open(index_file) as index_handle:
        index_handle.next()
        num_lines = int(index_handle.next().strip())
    print "Loading %d variants." % (num_lines)
    return num_lines

def grabix_index(fname):
    if not which("grabix"):
        print_cmd_not_found_and_exit("grabix")
    index_file = fname + ".gbi"
    if file_exists(index_file):
        return index_file
    print "Indexing {0} with grabix.".format(fname)
    subprocess.check_call("grabix index {fname}".format(fname=fname), shell=True)
    return index_file

def bgzip(fname):

    if not which("bgzip"):
        print_cmd_not_found_and_exit("bgzip")

    if is_gz_file(fname):
        return fname

    vcf_time = os.path.getmtime(fname)
    bgzip_file = fname + ".gz"

    if not file_exists(bgzip_file) or \
       (file_exists(bgzip_file) and os.path.getmtime(bgzip_file) < vcf_time):
        print "Bgzipping {0} into {1}.".format(fname, fname + ".gz")
        subprocess.check_call("bgzip -c {fname} > \
                              {fname}.gz".format(fname=fname),
                              shell=True)
    elif file_exists(bgzip_file) and os.path.getmtime(bgzip_file) > vcf_time:
        print "Loading with existing bgzip ({0}) version of {1}.".format(fname + ".gz", fname)

    return bgzip_file


def is_gz_file(fname):
    _, ext = os.path.splitext(fname)
    if ext == ".gz":
        return True
    else:
        return False

def get_submit_command(args):
    return "{cmd}"


def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    return os.path.exists(fname) and os.path.getsize(fname) > 0

def which(program):
    """ returns the path to an executable or None if it can't be found
     http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
     """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def combine_dicts(d1, d2):
    return dict(d1.items() + d2.items())

def get_ipython_args(args):
    return (args.scheduler, args.queue, args.cores)

def print_cmd_not_found_and_exit(cmd):
    sys.stderr.write("Cannot find {cmd}, install it or put it in your "
                     "path.".format(cmd=cmd))
    exit(1)

def use_scheduler(args):
    return bool(args.scheduler)
