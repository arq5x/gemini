import os
import warnings
from collections import namedtuple

import GeminiQuery

import tool_de_novo_mutations as de_novo_tool
import tool_autosomal_recessive as recessive_tool
import tool_autosomal_dominant as dominant_tool

# based upon bottle example here:
# https://bitbucket.org/timtan/bottlepy-in-real-case

# -- determine where I launch python and config lib path
# base_dir = os.path.dirname(__file__)
# third_party_path = os.path.abspath(os.path.join(base_dir, 'third_party' ))
# sys.path.insert(0, third_party_path)

# -- common bottle importation
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from bottle import TEMPLATE_PATH, Bottle, run, static_file, debug, request
    from bottle import jinja2_view as view, jinja2_template as template

debug(True)


base_dir = os.path.dirname(__file__)
TEMPLATE_PATH.append(os.path.abspath(os.path.join(base_dir, 'views')))

# -- the instance app is important
app = Bottle()

# -- serve static files, files located in static
static_folder = 'static'
_static_folder = os.path.join(os.path.dirname(__file__), static_folder)

@app.route('/stats/region/:chrom',method='GET')
def stats_region(chrom):
    # Note: chrom is give as an argument

    # we then extract start and end using HTML GET    
    start = request.GET.get('start', '').strip()
    end = request.GET.get('end', '').strip()

    # construct a query
    query =  "SELECT start, end from variants"
    query += " WHERE chrom = '" + chrom + "'"
    query += " AND start >= " + start
    query += " AND end <= " + end

    # issue the query
    gq = GeminiQuery.GeminiQuery(database)
    gq._set_gemini_browser(True)
    gq.run(query)

    # return query results in JSON format
    return{'features': [dict(row) for row in gq]}


@app.route('/static/<filepath:path>')
def server_static(filepath):
    return static_file(filepath, root=_static_folder)
# -- end of static folder configuration

# -- index page routing


@app.route('/index')
@app.route('/')
def index():
    return template('index.j2')

@app.route('/query_json', method='GET')
def query_json():
    query = request.GET.get('query', '').strip()
    
    gq = GeminiQuery.GeminiQuery(database)
    gq._set_gemini_browser(True)
    gq.run(query)
    
    return {'gemini_results': [dict(row) for row in gq]}


@app.route('/query', method='GET')
def query():

    def _get_fields():
        query = request.GET.get('query', '').strip()
        gt_filter = request.GET.get('gt_filter', '').strip()
        use_header = request.GET.get('use_header')
        igv_links = request.GET.get('igv_links')
        return query, gt_filter, use_header, igv_links
    
    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

        (query, gt_filter, use_header, igv_links) = _get_fields()

        if use_header: use_header = True
        if igv_links: igv_links = True

        gq = GeminiQuery.GeminiQuery(database)
        gq._set_gemini_browser(True)
        gq.run(query, gt_filter)


        if len(query) == 0:
            return template('query.j2', dbfile=database)

        if igv_links and ('chrom' not in query.lower()
                          or 'start' not in query.lower()
                          or 'end' not in query.lower()):
            return template('query.j2', dbfile=database,
                            rows=gq,
                            igv_links=igv_links,
                            igv_links_error=True,
                            use_header=use_header,
                            gt_filter=gt_filter,
                            query=query)
        else:
            return template('query.j2', dbfile=database,
                            rows=gq,
                            igv_links=igv_links,
                            igv_links_error=False,
                            use_header=use_header,
                            gt_filter=gt_filter,
                            query=query)

    # user clicked the "save to file" button
    elif request.GET.get('save', '').strip():

        (query, gt_filter, use_header, igv_links) = _get_fields()

        gq = GeminiQuery.GeminiQuery(database)
        gq.run(query, gt_filter)

        if len(query) == 0:
            return template('query.j2', dbfile=database)

        # dump the results to a text file.  this will be
        # stored in /static and a link will be given to
        # the user.
        tmp_file = '/tmp.txt'
        tmp = open(_static_folder + tmp_file, 'w')
        for row in gq:
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

    Arguments = namedtuple('Arguments', ['db', 'min_sample_depth'], verbose=True)
    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

        min_sample_depth = str(request.GET.get('min-depth', '').strip())
        igv_links = request.GET.get('igv_links')

        args = Arguments(db=database, min_sample_depth=min_sample_depth)

        de_novo_factory = \
            GeminiInheritanceModelFactory(args, model="de_novo")
        de_novo_factory.get_candidates()

        if len(min_sample_depth) == 0:
            row_iter = \
                de_novo_tool.get_de_novo_candidates(gq.c)
        else:
            row_iter = \
                de_novo_tool.get_de_novo_candidates(gq.c, int(min_sample_depth))

        return template('de_novo.j2', dbfile=database,
                        rows=row_iter,
                        igv_links=igv_links)

    else:
        return template('de_novo.j2', dbfile=database)


@app.route('/auto_rec', method='GET')
def auto_rec():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

        c = connect_to_db(database)

        row_iter = \
            recessive_tool.get_auto_recessive_candidates(c)

        return template('auto_rec.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_rec.j2', dbfile=database)


@app.route('/auto_dom', method='GET')
def auto_dom():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

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
