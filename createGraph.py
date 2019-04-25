from igraph import *

### USAGE
# python createGraph.py 3

### FORMATTING THE SEQUENCE

# inputs
with open('data.txt', 'r') as sequencefile:
    inputstr = sequencefile.read().replace('\n', '') # input string
k = int(sys.argv[1]) # k-mer of the reads

# split the string into reads
# here, I iterate through the input string by selecting k lengths
reads = [inputstr[i:i+k] for i in range(0, len(inputstr), k)]
# list of (k-1)-mer tuples for each read (left/right)
# here, I iterate through each read and split into k-1-mers, setting left and right
splitReads = [(read[0:k-1], read[1:]) for read in reads]
# list of (k-1)-mers, without repeats
# here, I iterate through each pair in splitReads, and a python set prevents repeats
vertexLabels = list(set([read for pair in splitReads for read in pair]))

# vertexLabels = ["AA", "AB", "BB", "BA"]
# splitReads = [("AA", "AA"), ("AA", "AB"), ("AB", "BB"), ("BB", "BB"), ("BB", "BA")]

### TESTING (DELETE)
print(vertexLabels)
print(splitReads)

### CREATING THE GRAPH

# g = Graph()
# g.add_vertices(len(vertexLabels))
# g.vs["seq"] = vertexLabels
# for pair in splitReads:
#     left = g.vs.select(seq_eq = pair[0])     #not sure about the runtime of this - might have to store a deparate map from vertexLabel -> vertexID
#     right = g.vs.select(seq_eq = pair[1])
#     g.add_edges(c(left, right))

g = Graph()
g.add_vertices(len(vertexLabels))
g.add_edges([1,2])
print(g)
