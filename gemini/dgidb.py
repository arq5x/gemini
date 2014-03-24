import urllib2
import json

def query_dgidb(genes):
    """
    Batch query DGIdb for drug-gene interaction data for
    a set of genes.
    """

    def convert(input):
        """
        Convert JSON UNICODE to plain ole strings.
        """
        if isinstance(input, dict):
            return {convert(key): convert(value) for key, value in input.iteritems()}
        elif isinstance(input, list):
            return [convert(element) for element in input]
        elif isinstance(input, unicode):
            return input.encode('utf-8')
        else:
            return input

    # make a single request to DGIdb for all of the genes requested
    dgidb_url = 'http://dgidb.genome.wustl.edu/api/v1/interactions.json?genes='
    
    # None is present by default. Make sure we have more than None
    if len(genes) > 1:
        query = dgidb_url + ','.join(genes.keys())
        
        response = urllib2.urlopen(query)
        data = convert(json.load(response))
        matches = data['matchedTerms']

        # store the results for all of the genes. if there are no matches
        # in DGIdb, the result will be None.
        gene_dgidb_info = {}
        for gene in genes:
            gene_dgidb_info[gene] = None

        for match in matches:
            gene = match['searchTerm']
            gene_dgidb_info[gene] = dict(match)

        return gene_dgidb_info
    else:
        return None