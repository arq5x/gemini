#!/usr/bin/env python
import sqlite3
import os
import database as gemini_db
import gemini_utils as util

def merge_db_chunks(args):


    # open up a new database
    if os.path.exists(args.db):
        os.remove(args.db)
    
    main_conn = sqlite3.connect(args.db)
    main_conn.isolation_level = None
    main_conn.row_factory = sqlite3.Row
    main_curr = main_conn.cursor()
    main_curr.execute('PRAGMA synchronous = OFF')
    main_curr.execute('PRAGMA journal_mode=MEMORY')
    # create the gemini database tables for the new DB
    gemini_db.create_tables(main_curr)
    
    
    databases = []
    for database in args.chunkdbs:
        databases.append(database)

    for idx, database in enumerate(databases):
    
        cmd = "attach ? as toMerge"
        main_curr.execute(cmd, (database[0], ))
        
        cmd = "INSERT INTO variants select * FROM toMerge.variants"
        main_curr.execute(cmd)
        
        cmd = \
        "INSERT INTO variant_impacts select * FROM toMerge.variant_impacts"
        main_curr.execute(cmd)
        
        # we only need to add these tables from one of the chunks.
        if idx == 0:
            cmd = "INSERT INTO samples select * FROM toMerge.samples"
            main_curr.execute(cmd)

            cmd = "INSERT INTO resources select * FROM toMerge.resources"
            main_curr.execute(cmd)

            
            
        cmd = "detach toMerge"
        main_curr.execute(cmd)

    gemini_db.create_indices(main_curr)
    main_conn.commit()
    main_curr.close

def merge_chunks(parser, args):
    merge_db_chunks(args)
