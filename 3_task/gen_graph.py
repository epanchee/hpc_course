from __future__ import division

import networkx as nx
import numpy as np

G = nx.erdos_renyi_graph(30, 0.05)
# G = nx.Graph()
# G.add_edges_from([(1,0),(2,0),(3,0),(4,0),(1,5),(5,4)])
M = nx.to_scipy_sparse_matrix(G)

f = open('csr_graph', 'w')

for arr in [M.data, M.indices, M.indptr]:
	# # f.write("%d\n%s\n" % (len(arr), " ".join([str(x-1) if x !=0 else str(x) for x in arr])))
	f.write("%d\n%s\n" % (len(arr), " ".join([str(x) for x in arr])))

# f.write("%d\n%s\n" % (len(M.data), " ".join([str(x) for x in M.data])))
# f.write("%d\n%s\n" % (len(M.indices), " ".join([str(x-1) if x !=0 else str(x) for x in M.indices])))
# f.write("%d\n%s\n" % (len(M.indptr), " ".join([str(x) for x in M.indptr])))
	
# print "data %r indices %r indptr %r" % (M.data, M.indices, M.indptr)

# checking

V = G.nodes()
P = [[] for _ in xrange(G.number_of_nodes())]
BC = [0]*G.number_of_nodes()
sigma = [0]*G.number_of_nodes()
d = [0]*G.number_of_nodes()
delta = [0]*G.number_of_nodes()

for s in V:
	S = [] # stack
	P = [[] for _ in xrange(G.number_of_nodes())]
	sigma = [0]*G.number_of_nodes()
	d = [-1]*G.number_of_nodes()
	sigma[s] = 1
	d[s] = 0
	Q = []
	Q.append(s) # push to queue
	
	while Q:
		v = Q.pop(0)
		S.append(v)
		for w in G.neighbors(v):
			if d[w] < 0:
				Q.append(w)
				d[w] = d[v] + 1
			if d[w] == d[v] + 1:
				# print "s %d sigma[w] %d sigma[v] %d" % (s, sigma[w], sigma[v])
				sigma[w] = sigma[w] + sigma[v]
				P[w].append(v)
	
	delta = [0]*G.number_of_nodes()
	
	while S:
		w = S.pop()
		# print "w: %d size of P[w]: %d" % (w, len(P[w]))
		for v in P[w]:
			# print "v: %d delta[v]: %d sigma[v]: %d sigma[w]: %d sigma[v]/sigma[w]: %.3f" % (v, delta[v], sigma[v], sigma[w], sigma[v]/sigma[w])
			delta[v] = delta[v] + sigma[v]/sigma[w] * (1 + delta[w])
		if w != s:
			BC[w] = BC[w] + delta[w]
	
f = open('btwcheck_py', 'w')
	
BC = ["%.2f" % (x/2) for x in BC]
print BC
print "NetworkX algorithm:"
print nx.betweenness_centrality(G, normalized=False)

f.write(" ".join(BC))
f.write("\n")
f.write(str(nx.betweenness_centrality(G, normalized=False)))