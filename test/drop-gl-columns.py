"""
This is used for testing backwards compat.
It drops the genotype likelihood columns from the db.
"""

import sqlite3
import sys
import shutil

dbfile = sys.argv[1]
newname = sys.argv[2]
shutil.copyfile(dbfile, newname)

db = sqlite3.connect(newname)

cur = db.execute("PRAGMA table_info('variants')")

# name, type
cols = [(row[1], row[2]) for row in cur if not "phred" in row[1]]

db.execute("BEGIN")
db.execute("CREATE TABLE v2(%s)" % (",\n".join("%s %s" % c for c in cols)))
db.execute("INSERT INTO  v2 select %s from variants" % (",".join(c[0] for c in cols)))
db.execute("DROP TABLE variants")
db.execute("CREATE TABLE variants(%s)" % (",\n".join("%s %s" % c for c in cols)))
db.execute("INSERT INTO  variants select * from v2")
db.execute("DROP TABLE v2")
