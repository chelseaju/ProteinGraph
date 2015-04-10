from echo import *
"""
	get edge label from reference
"""
def get_edge_label(dist, ref):
	for (e, l) in ref:
		if(float(dist) <= e):
			return l

	return None
"""
	Build the edge label guid
"""

def build_edge_guideline(edge_file):

	echo("Building edge guideline using %s" %(edge_file))

	edge_labels = []
	fh = open(edge_file, 'r')
	for line in fh:
		(threshold, label) = line.rstrip().split('\t')
		edge_labels.append((float(threshold), label))
	fh.close()

	return edge_labels


"""
	convert pdb to graph
"""
def pdb_to_graph(dist_info, graph_file, edge_labels, graph_id, graph_title):

	edges = []
	vertices = {}

	e_count = 1 
	for (a, b, dist) in dist_info:
		e_label = get_edge_label(dist, edge_labels)
		
		(a_label, a_id) = a.split("_")
		(b_label, b_id) = b.split("_")

		vertices[int(a_id)] = a_label
		vertices[int(b_id)] = b_label
	
		if(e_label):
			edges.append(( int(a_id), int(b_id), e_label, e_count))
			e_count = e_count + 1

	#write to file
	fh = open(graph_file, 'w')
	fh.write("t # %d %d %s\n" %(graph_id, len(vertices.keys()), graph_title))
	for key in sorted(vertices.keys()):
		fh.write("v %d %s\n" %(key, vertices[key]))


	for e in edges:
		fh.write("e %d %d %s %d\n" %(e[0], e[1], e[2], e[3]))

	fh.close()

	echo("Writing graph to %s" %(graph_file))

"""
	retrieve pdb pairwise distance from file
"""
def pdb_file_to_graph(dist_file, graph_file, edge_labels):

	echo("Parsing distance infor from %s" %(dist_file))

	dist_info = []
	fh = open(dist_file, 'r')
	for line in fh:
		(a, b, dist) = line.rstrip().split()
		dist_info.append((a, b, dist))
	fh.close()

	pdb_to_graph(dist_info, graph_file, edge_labels, 0, "title")




if __name__ == "__main__":
	edge_ref = "edge_info_v2.txt"
	pdb_dist = "test.txt"
	graph_file = "graph.txt"

	edge_labels = build_edge_guideline(edge_ref)
	pdb_file_to_graph(pdb_dist, graph_file, edge_labels)

