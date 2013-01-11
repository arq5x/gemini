import sys
import os
import logging
import sqlite3
import copy
from gemini_query import filter_query, \
                         apply_basic_query, \
                         apply_query_w_genotype_select

import tool_de_novo_mutations as de_novo_tool
import tool_autosomal_recessive as recessive_tool
import tool_autosomal_dominant as dominant_tool

# based upon bottle example here:
# https://bitbucket.org/timtan/bottlepy-in-real-case

# -- determine where I launch python and config lib path
#base_dir = os.path.dirname(__file__)
#third_party_path = os.path.abspath(os.path.join(base_dir, 'third_party' ))
#sys.path.insert(0, third_party_path)

# -- common bottle importation
from bottle import TEMPLATE_PATH, Bottle, run, CGIServer, static_file, debug, request
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


def connect_to_db(database):

    conn = sqlite3.connect(database)
    conn.isolation_level = None
    conn.row_factory = sqlite3.Row # allow us to refer to columns by name
    c = conn.cursor()
    
    return c

# -- query routing
def process_query(query, gt_filter, use_header):
    
    c = connect_to_db(database)
        
    if query and gt_filter:
        row_iter = filter_query(c, query, gt_filter, use_header)
    else:
        query_pieces = query.split()
        if not any(s.startswith("gt") for s in query_pieces) and \
           not any("gt" in s for s in query_pieces):
            row_iter = apply_basic_query(c, query, use_header)
        else:
            row_iter = apply_query_w_genotype_select(c, query, use_header)
    
    return row_iter

@app.route('/query', method='GET')
def query():
    
    def _get_fields():
        query      = request.GET.get('query', '').strip()
        gt_filter  = request.GET.get('gt_filter', '').strip()
        use_header = request.GET.get('use_header')
        igv_links  = request.GET.get('igv_links')
        
        row_iter = process_query(query, gt_filter, use_header)
        
        print use_header 
        return query, gt_filter, use_header, igv_links, row_iter

    # user clicked the "submit" button
    if request.GET.get('submit','').strip():

        (query, gt_filter, use_header, igv_links, row_iter) = _get_fields()
        
        if len(query) == 0: return template('query.j2', dbfile=database)
        
        if igv_links and ('chrom' not in query.lower() \
                          or 'start' not in query.lower() \
                          or 'end' not in query.lower()):
            return template('query.j2', dbfile=database, 
                            rows=row_iter,
                            igv_links=igv_links,
                            igv_links_error=True,
                            use_header=use_header,
                            gt_filter=gt_filter,
                            query=query)
        else:
            return template('query.j2', dbfile=database, 
                            rows=row_iter,
                            igv_links=igv_links,
                            igv_links_error=False,
                            use_header=use_header,
                            gt_filter=gt_filter,
                            query=query)
                                    
    # user clicked the "save to file" button
    elif request.GET.get('save','').strip():
        
        (query, gt_filter, use_header, igv_links, row_iter) = _get_fields()
        
        if len(query) == 0: return template('query.j2', dbfile=database)

        # dump the results to a text file.  this will be
        # stored in /static and a link will be given to 
        # the user.
        tmp_file =  '/tmp.txt'
        tmp = open(_static_folder + tmp_file, 'w')
        for row in row_iter:
            tmp.write('\t'.join(str(c) for c in row) + '\n')
        tmp.close()
        
        return template('query.j2', dbfile=database, 
                                    tmp_file=tmp_file,
                                    igv_links=igv_links,
                                    igv_links_error=True,
                                    use_header=use_header,
                                    gt_filter=gt_filter,
                                    query=query)
    # user did nothing.
    else:
        return template('query.j2', dbfile=database)


@app.route('/de_novo', method='GET')
def de_novo():
    
    # user clicked the "submit" button
    if request.GET.get('submit','').strip():
        
        min_sample_depth  = str(request.GET.get('min-depth', '').strip())
        igv_links = request.GET.get('igv_links')

        c = connect_to_db(database)
    
        if len(min_sample_depth) == 0:
            row_iter = \
                de_novo_tool.get_de_novo_candidates(c)
        else: 
            row_iter = \
                de_novo_tool.get_de_novo_candidates(c, int(min_sample_depth))

        return template('de_novo.j2', dbfile=database, 
                                      rows=row_iter, 
                                      igv_links=igv_links)

    else:
        return template('de_novo.j2', dbfile=database)
        
@app.route('/auto_rec', method='GET')
def auto_rec():
    
    # user clicked the "submit" button
    if request.GET.get('submit','').strip():

        c = connect_to_db(database)

        row_iter = \
            recessive_tool.get_auto_recessive_candidates(c)

        return template('auto_rec.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_rec.j2', dbfile=database)


@app.route('/auto_dom', method='GET')
def auto_dom():
    
    # user clicked the "submit" button
    if request.GET.get('submit','').strip():

        c = connect_to_db(database)

        row_iter = \
            dominant_tool.get_auto_dominant_candidates(c)

        return template('auto_dom.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_dom.j2', dbfile=database)


@app.route('/db_schema', method='GET')
def db_schema():
    return template('db_schema.j2')

def browser_main(parser, args):
    global database 
    
    print "!!!!!"
    print "NOTE: open a browser and point it to http://localhost:8088/query"
    print "!!!!!"
    
    database = args.db
    
    run(app, host='localhost', port=8088, 
             reloader=True, debug=True)


