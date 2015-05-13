#!/usr/bin/env python

############################################################
# Gives a final pathway map of uniprot -> hsa -> gene for ensembl
#########################################################
from collections import defaultdict

filename = 'map_uni2hsa2gene'
out = open(filename, 'w')

filep = 'add_path'
outf = open(filep, 'w') 

unihsa = defaultdict(list)# dict of uniprot to hsa map
hsagene = defaultdict(list) # dict of hsa to gene map
map = defaultdict(list) # 
pathway = defaultdict(list)
add = defaultdict(list)


for line in open("uniprot_hsa_map", 'r'):
    r = line.strip().split("\t")
    (key, value) = (r[0], r[1])
    unihsa[key].append(value)

for each in open("kegg_hsa_to_gene_fmt", 'r'):
    l = each.strip().split("\t")
    (key, value) = (l[0], l[1])
    hsagene[key] = value

prim_value_is_key = [(k, each, hsagene[each]) for k, v in unihsa.iteritems() for each in v if each in hsagene]
for each in prim_value_is_key:
    genes = each[2].split(",")
    for gene in genes:
        g = gene.lstrip() # remove leading white spaces before gene names
        string = [str(each[0]), str(each[1]), str(g)]
        out.write("\t".join(string))
        out.write('\n')
out.close()

for line in open("map_uni2hsa2gene", 'r'):
    p = line.strip().split("\t")
    (key, value) = ((str(p[0]), str(p[2])), str(p[1])) 
    map[key].append(value)

for path in open("kegg_pathway", 'r'):
    k = path.strip().split("\t")
    (key, value) = (k[0], k[1])
    pathway[key] = value

for lines in open("kegg_hsa_to_path", 'r'):
    u = lines.strip().split("\t")
    hs = str(u[0])
    pathcode = str(u[1])
    thread = [hs, pathcode, pathway[pathcode]]
    outf.write("\t".join(thread))
    outf.write("\n")
outf.close()

for rows in open("add_path", 'r'):
    t = rows.strip().split("\t")
    (key, value) = (str(t[0]), (str(t[1]), str(t[2])))
    add[key].append(value)
    
for line in open("ensembl_genes66", 'r'):
    k = line.strip().split("\t")
    uniprot = str(k[0])
    agn = str(k[1])
    hgnc = str(k[2])
    ens_geneid = str(k[3])
    ens_transid = str(k[4])
    
    hsa = [(uniprot, agn, hgnc, ens_geneid, ens_transid, v) for k, v in map.iteritems() if (k[1] == agn and k[0] == uniprot) or (k[1] == hgnc and k[0] == uniprot)]
    if not hsa:
        print uniprot + "\t" + agn + "\t" + hgnc + "\t" + ens_geneid + "\t" + ens_transid + "\t" + "None" + "\t" + "None"   
    
    elif hsa: # if list is not empty
        for each in hsa:
            for keggid in each[5]:
                if keggid in add:
                    for eachvalue in add[keggid]:
                        pathstring = ";".join(eachvalue)
                        print str(each[0]) + "\t" + str(each[1]) + "\t" + str(each[2]) + "\t" + str(each[3]) + "\t" + str(each[4]) + "\t" + str(keggid) + "\t" + pathstring
                else:
                    print str(each[0]) + "\t" + str(each[1]) + "\t" + str(each[2]) + "\t" + str(each[3]) + "\t" + str(each[4]) + "\t" + str(keggid) + "\tNone"
