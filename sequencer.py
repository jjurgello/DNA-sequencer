# Jennifer Steele (jsteele3), JJ Urgello (urgello2)
# CS466 Final Project, DNA Sequencer

import igraph
import sys
import re
from time import time

### USAGE
# python sequencer.py [k, not necessary]

def piggysequency(k, n):
    
    ### FORMATTING THE SEQUENCE

    # input formating
    validchars = ['A', 'T', 'C', 'G']
    with open('fullData.txt', 'r') as sequencefile: # set inputstr
        inputstr = ''
        while len(inputstr) < n:
            c = sequencefile.read(1)
            if c in validchars:
                inputstr += c
    
    # split the string into reads
    # here, I iterate through the input string by selecting k lengths
    reads = [inputstr[i:i+k] for i in range(0, len(inputstr) - (k-1))]
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

    # display graph
    # print(g)
    print('n: ', len(vertexLabels))
    print('m: ', len(edges))
    
    ### SEQUENCING
    
    # Hierholzer's Algorithm to implement to construction of a Eulerian Path
    # reference: https://iampandiyan.blogspot.com/2013/10/c-program-to-find-euler-path-or-euler.html
    
    numEdges = len(splitReads) # count edges
    currVertex = g.vs.find(splitReads[0][0])   # OR find starting vertex - odd outward degree
    breadCrumbs = [currVertex] # temporary path, order of verticies visited
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
    
    lastVertex = breadCrumbs.pop()
    yellowBrickRd.append(lastVertex['name'])
    yellowBrickRd.reverse()
    
    ### RESULT
    
    # construct the sequence
    result = yellowBrickRd[0]
    result = result[:-1]
    for vertexStr in yellowBrickRd:
        result += vertexStr[-1:]
    
    return inputstr, result

### RUNNING and ANALYSIS

# a default run
k = int(sys.argv[1]) if len(sys.argv) > 1 else 3 # k-mer of the reads
n = int(sys.argv[2]) if len(sys.argv) > 2 else 20 # length of sequence to read 
ins, res = piggysequency(k, n)
print("input:  ", ins)
print("result: ", res)

# for analysis
if len(sys.argv) <= 1:
    # comparing different length sequences
    t0 = time()
    piggysequency(3, 10) # 10 letters
    t1 = time()
    piggysequency(3, 100) # 100 letters
    t2 = time()
    piggysequency(3, 1000) # 1000 letters
    t3 = time()
    piggysequency(3, 10000) # 10000 letters
    t4 = time()
    piggysequency(3, 100000) # 100000 letters
    t5 = time()

    print('10 letter runtime: ', t1-t0)
    print('100 letter runtime: ', t2-t1)
    print('1000 letter runtime: ', t3-t2)
    print('10000 letter runtime: ', t4-t3)
    print('100000 letter runtime: ', t5-t4)

    # comparing different k-mer lengths
    t0 = time()
    piggysequency(3, 1000) # k = 3
    t1 = time()
    piggysequency(9, 1000) # k = 3
    t2 = time()
    piggysequency(27, 1000) # k = 3
    t3 = time()
    piggysequency(81, 1000) # k = 3
    t4 = time()
    piggysequency(243, 1000) # k = 3
    t5 = time()
    
    print('k = 3 runtime: ', t1-t0)
    print('k = 9 runtime: ', t2-t1)
    print('k = 27 runtime: ', t3-t2)
    print('k = 81 runtime: ', t4-t3)
    print('k = 243 runtime: ', t5-t4)


