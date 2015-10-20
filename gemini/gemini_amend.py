from gemini_subjects import get_subjects
from ped import load_ped_file, get_ped_fields
from gemini_utils import quote_string
from database import database_transaction

def amend(parser, args):
    if args.db is None:
        parser.print_help()
        exit("ERROR: amend needs a database file.")
    if args.sample:
        amend_sample(args)

def amend_sample(args):
    loaded_subjects = get_subjects(args)
    ped_dict = load_ped_file(args.sample)
    header = get_ped_fields(args.sample)
    with database_transaction(args.db) as c:
        add_columns(header, c, args.clear)
        for k, v in loaded_subjects.items():
            if k in ped_dict:
                item_list = map(quote_string, ped_dict[k])
                sample = zip(header, item_list)
                set_str = ",".join([str(x) + "=" + str(y) for (x, y) in sample])
                sql_query = "update samples set {0} where sample_id={1}"
                c.execute(sql_query.format(set_str, v.sample_id))

def add_columns(header, c, clear=False):
    """
    add any missing columns to the samples table
    """
    for column in header:
        try:
            c.execute('ALTER TABLE samples ADD COLUMN {0}'.format(column))
        except:
            pass
        if clear:
            c.execute("UPDATE samples set {0} = NULL;".format(column))
