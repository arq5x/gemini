import os
import warnings
import webbrowser
from collections import namedtuple

import GeminiQuery

import GeminiQuery
from gemini.gim import (AutoDom, AutoRec, DeNovo, MendelViolations, CompoundHet)

database = None

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
    from bottle import jinja2_template as template

debug(True)


base_dir = os.path.dirname(__file__)
TEMPLATE_PATH.append(os.path.abspath(os.path.join(base_dir, 'views')))

# -- the instance app is important
app = Bottle()

# -- serve static files, files located in static
static_folder = 'static'
_static_folder = os.path.join(os.path.dirname(__file__), static_folder)

@app.route('/stats/region/:chrom', method='GET')
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


default_cols = ['chrom', 'start', 'end', 'ref', 'alt',
                'polyphen_pred', 'sift_pred',
                'max_aaf_all', 'impact', 'impact_severity',
                'gene', 'biotype', 'transcript']
# turn a dictionary into something that can be accessed by attribute.
class Arguments(object):
    """
    >>> args = Arguments(db='some.db')
    >>> args.db
    'some.db'
    """
    _opts = ("columns", "db", "filter", "min_kindreds", "families",
                 "pattern_only", "max_priority", # only for compound_het
                 "allow_unaffected", "min_gq", "lenient", "min_sample_depth")
    def __init__(self, **kwargs):
        if not "min_gq" in kwargs: kwargs['min_gq'] = 0
        if not "lenient" in kwargs: kwargs['lenient'] = False
        for k in ("families", "filter"):
            if not k in kwargs: kwargs[k] = None
        if not "gt_phred_ll" in kwargs: kwargs['gt_phred_ll'] = None
        if not "min_sample_depth" in kwargs: kwargs['min_sample_depth'] = 0

        for k in ("min_kindreds", "max_priority"):
            if not k in kwargs: kwargs[k] = 1
        for k in ("pattern_only", "allow_unaffected"):
            if not k in kwargs: kwargs[k] = False
        if not "columns" in kwargs:
            kwargs['columns'] = ",".join(default_cols)
        self.__dict__.update(**kwargs)


@app.route('/de_novo', method='GET')
def de_novo():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():

        min_sample_depth = str(request.GET.get('min-depth', '10').strip())
        igv_links = request.GET.get('igv_links')

        args = Arguments(db=database, min_sample_depth=min_sample_depth)

        row_iter = DeNovo(args).report_candidates()

        return template('de_novo.j2', dbfile=database,
                        rows=row_iter,
                        igv_links=igv_links)

    else:
        return template('de_novo.j2', dbfile=database)


@app.route('/auto_rec', method='GET')
def auto_rec():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():
        min_sample_depth = str(request.GET.get('min-depth', '10').strip())
        args = Arguments(db=database, min_sample_depth=min_sample_depth)

        row_iter = AutoRec(args).report_candidates()
        return template('auto_rec.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_rec.j2', dbfile=database)


@app.route('/auto_dom', method='GET')
def auto_dom():

    # user clicked the "submit" button
    if request.GET.get('submit', '').strip():
        min_sample_depth = str(request.GET.get('min-depth', '10').strip())
        args = Arguments(db=database, min_sample_depth=min_sample_depth)
        row_iter = AutoDom(args).report_candidates()
        return template('auto_dom.j2', dbfile=database, rows=row_iter)

    else:
        return template('auto_dom.j2', dbfile=database)


@app.route('/db_schema', method='GET')
def db_schema():
    return template('db_schema.j2')


## Switch between the different available browsers
def browser_puzzle(args):
    host = args.host
    port = args.port

    try:
        # Puzzle browser plugin
        from puzzle.server import factory as puzzle_app
        from puzzle.plugins import GeminiPlugin
        from puzzle.server.settings import BaseConfig
    except ImportError as e:
        raise ImportError("%s\nPlease 'pip install puzzle' if you want to run it\n" % e)

    plugin = GeminiPlugin(db=args.db, vtype="sv")
    root = os.path.expanduser("~/.puzzle")

    BaseConfig.PUZZLE_BACKEND = plugin
    BaseConfig.UPLOAD_DIR = os.path.join(root, 'resources')

    puzzle_srv = puzzle_app.create_app(config_obj=BaseConfig)
    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(puzzle_srv, host=host, port=port)

def browser_builtin(args):
    host = args.host
    port = args.port

    webbrowser.open_new_tab("http://{}:{}".format(host, port))
    run(app, host=host, port=port,
        reloader=True, debug=True)


def browser_main(parser, args):
    global database
    database = args.db
    browser = args.use

    try:
        if args.use == "puzzle":
            browser_puzzle(args)
        # XXX: https://github.com/dgaston/kvasir
        #if args.use == "kvasir":
        #    raise NotImplementedError
        elif args.use == "builtin":
            browser_builtin(args)
        else:
            raise NotImplementedError("GEMINI-compatible Browser '{browser}' not found.".format(browser=browser))
    except ImportError as e:
        raise ImportError("{exc}\nIs {browser} correctly installed?".format(exc=e, browser=browser))
