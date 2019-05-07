# jennifer steele, jj urgello
# CS466 final project, DNA sequencer

import igraph
import sys
import re

### USAGE
# python sequencer.py [k, not necessary]

### FORMATTING THE SEQUENCE

# inputs
validchars = ['A', 'T', 'C', 'G']
with open('data.txt', 'r') as sequencefile: # set inputstr
    inputstr = ''.join([c for c in sequencefile.read() if c in validchars])
k = int(sys.argv[1]) if len(sys.argv) > 1 else 3 # k-mer of the reads

# split the string into reads
# here, I iterate through the input string by selecting k lengths
reads = [inputstr[i:i+k] for i in range(0, len(inputstr) - k)]
#if len(inputstr) % k != 0:
#    reads = reads[0:len(reads)-1]
# list of (k-1)-mer tuples for each read (left/right)
# here, I iterate through each read and split into k-1-mers, setting left and right
splitReads = [(read[0:k-1], read[1:]) for read in reads]
# list of (k-1)-mers, without repeats
# here, I iterate through each pair in splitReads, and a python set prevents repeats
vertexLabels = list(set([read for pair in splitReads for read in pair]))

### CREATING THE GRAPH

# create an igraph
g = igraph.Graph()
# add the vertex labels
g.add_vertices(vertexLabels)
# add the edges, determined by splitReads
edges = []
for l, r in splitReads:
    edges.append((g.vs.find(l), g.vs.find(r)))
g.add_edges(edges)

# display the graph
print(g)

### SEQUENCING

#
