import sys
import os
import logging
import sqlite3
from gemini_query import filter_query, \
                         apply_basic_query, \
                         apply_query_w_genotype_select

# based upon bottle example here:
# https://bitbucket.org/timtan/bottlepy-in-real-case

# -- determine where I launch python and config lib path
#base_dir = os.path.dirname(__file__)
#third_party_path = os.path.abspath(os.path.join(base_dir, 'third_party' ))
#sys.path.insert(0, third_party_path)

# -- common bottle importation
from bottle import TEMPLATE_PATH, Bottle, run, static_file, debug, request
from bottle import jinja2_view as view, jinja2_template as template
debug(True)


base_dir = os.path.dirname(__file__)
TEMPLATE_PATH.append(os.path.abspath(os.path.join(base_dir, 'views' )))

# -- the instance app is important
app = Bottle()

# -- serve static files, files located in static 
static_folder = 'static'
_static_folder = os.path.join( os.path.dirname(__file__), static_folder)
@app.route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root=_static_folder)
# -- end of static folder configuration

# -- index page routing
@app.route('/index')
@app.route('/')
def index():
    return template('index.j2')

# -- query routing
@app.route('/query', method='GET')
def query():

    if request.GET.get('submit','').strip():

        query     = request.GET.get('query', '').strip()
        gt_filter = request.GET.get('gt-filter', '').strip()
        header    = request.GET.get('header')

        if len(query) == 0:
            return template('query.j2')

        conn = sqlite3.connect(database)
        conn.isolation_level = None
        conn.row_factory = sqlite3.Row # allow us to refer to columns by name
        c = conn.cursor()
        
        if query and gt_filter:
            row_iter = filter_query(c, query, gt_filter, header)
        else:
            query_pieces = query.split()
            if not any(s.startswith("gt") for s in query_pieces) and \
               not any("gt" in s for s in query_pieces):
                row_iter = apply_basic_query(c, query, header)
            else:
                row_iter = apply_query_w_genotype_select(c, query, header)

        return template('query.j2', rows=row_iter, query=query)

    else:
        return template('query.j2')


def browser_main(parser, args):
    global database 
    
    print "!!!!!"
    print "NOTE: open a browser and point it to http://localhost:8088/query"
    print "!!!!!"
    
    database = args.db
    run(app, host='localhost', port=8088, reloader=True, debug=True)



