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
g.to_directed()
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

# Hierholzer's Algorithm to implement to construction of a Eulerian Path
# reference: https://iampandiyan.blogspot.com/2013/10/c-program-to-find-euler-path-or-euler.html

numEdges = len(splitReads) # count edges
currVertex = g.vs.find(vertexLabels[0])   # OR find starting vertex - odd outward degree
breadCrumbs = [] # temporary path, order of verticies visited
yellowBrickRd = [] # backwards final path

while len(yellowBrickRd) < numEdges:
    # set currEdge to random non-visited outgoing edge from currVertex
    neighs = g.neighbors(currVertex, mode=1)

    if len(neighs) is 0: # if there are no more edges to visit
        # set currVertex to pop from breadCrumbs
        currVertex = breadCrumbs.pop()
        # add currVertex to yellowBrickRd
        yellowBrickRd.append(g.vs[currVertex]['name'])
    else:
        # pick a random edge from that list
        currEdge = g.get_eid(currVertex, neighs[0])
        currVertex = neighs[0]
        # add currVertex to breadCrumbs
        breadCrumbs.append(currVertex)
        # delete edge from graph
        g.delete_edges(currEdge)

yellowBrickRd.reverse()

### RESULT

# print the sequence list
print(yellowBrickRd)

# construct the sequence
#result = ''
#for s in yellowBrickRd:
#J    ''.join(s

print(inputstr)
