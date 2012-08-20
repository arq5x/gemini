#!/usr/bin/env python

import cPickle
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
from pygraph.classes.exceptions import AdditionError
gr = graph()

output = open('hprd_interaction_graph', 'wb')
for line in open("BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", 'r'):
    fields=line.strip().split("\t")
    first = str(fields[0])
    second = str(fields[3])
    if first != "-":
        try:
            gr.add_nodes([first])
        except AdditionError:
            pass;
    if second != "-":
        try:
            gr.add_nodes([second])
        except AdditionError:
            pass;

    if (first == "-" or second == "-" or first == second):
        pass;
    else:
        try:
            gr.add_edge((first, second))
        except AdditionError:
            pass;

cPickle.dump(gr, output)
output.close()
              
