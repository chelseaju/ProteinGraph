"""
	Select positive and negative sets of protein for a given family
"""


import argparse, re, os, random
from pdb import *
from graph import *
from echo import *
from datetime import datetime

# negative graphs
def retrieve_negative(pdb_file, edge_info, num, dir):

	candidates = []
	fh = open(pdb_file, 'r')
	for line in fh:
		pdb_id = line.rstrip().split()[0]
		candidates.append(pdb_id)
	fh.close()

	random.seed(datetime.now())
	selected_id = random.sample(candidates, int(num))

	negative_count = 0
	for pdb_chain in sorted(selected_id):
		pdb_id = pdb_chain[:-1]
		chain = pdb_chain[-1]

		# distance file name
		pdb_dst_file = dir + pdb_id + "_" + chain + ".dist"
		
		# graph file name
		pdb_graph_file = dir + pdb_id + "_" + chain + ".txt"

		# parse pdb data
		pdb_info = parse_pdb_by_id(pdb_id, chain)
	
                if(pdb_info):
                        # comput distance
                        pw_dist = pairwise_distance(pdb_info, pdb_dst_file)

                        # convert structure to graph
			title = pdb_id + "_" + chain
                        pdb_to_graph(pw_dist, pdb_graph_file, edge_info, negative_count, title)

                        negative_count = negative_count + 1

                if(negative_count % 100 == 0):
                        time.sleep(3)

	return negative_count
		

def main(parser):
	options = parser.parse_args()
	edge_ref = options.eref
	pdb = options.pdb	
	num = options.num
	dir = options.dir	


	if(dir[-1] != "/"):
		dir += "/"

	os.system("mkdir -p %s" %(dir))

	# build references
	edge_info = build_edge_guideline(edge_ref)

	negative_count = retrieve_negative(pdb, edge_info, num, dir)

	# merge negative and positive graph
	os.system("echo %d > %stmp.txt" %(negative_count, dir))
	os.system("cat %s*_*.txt >> %stmp.txt" %(dir, dir))

        aa_to_number( "%stmp.txt" %(dir), "%snegative_graphs.txt" %(dir), [2])

	echo("Writing negative graph to %snegative_graphs.txt" %(dir))

	echo("Done")

if __name__ == "__main__":

        parser = argparse.ArgumentParser(prog='03_select_negative_sets.py')
        parser.add_argument("-e", "--edge_reference", dest = "eref", type=str, help="edge reference file", required = True)
	parser.add_argument("-p", "--pdb", dest= "pdb", type=str, help="filtered pdb", required = True)
	parser.add_argument("-n", "--num", dest= "num", type=str, help="number of proteins", required = True)
	parser.add_argument("-d", "--directory", dest = "dir", type=str, help="directory for output", required = True)
	main(parser)


