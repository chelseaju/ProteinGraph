"""
	Select positive set of proteins with a given family id
"""


import argparse, re, os, time, random
from pdb import *
from graph import *
from echo import *

# positive graphs
def retrieve_positive(fam_ref, fam_id, edge_info, dir):
	positive_count = 0
	POSITIVE = []

	fh = open(fam_ref, 'r')
	for line in fh:
		(pdb_id, family) = line.rstrip().split()
		if(family == fam_id):
			POSITIVE.append(pdb_id)
	fh.close()

	POSITIVE = set(POSITIVE)

	for pdb_chain in sorted(POSITIVE):
		pdb = pdb_chain[:-1]
		chain = pdb_chain[-1]

		# distance file name
		pdb_dst_file = dir + pdb + "_" + chain + ".dist"

		# graph file name
		pdb_graph_file = dir + pdb + "_" + chain + ".txt"

		# parse pdb data
		pdb_info = parse_pdb_by_id(pdb, chain)

		if(pdb_info):
			# comput distance
			pw_dist = pairwise_distance(pdb_info, pdb_dst_file)
	
			# convert structure to graph
			title = fam_id+ " " + pdb_chain
			pdb_to_graph(pw_dist, pdb_graph_file, edge_info, positive_count, title)
		
			positive_count = positive_count + 1

		if(positive_count % 100 == 0):
			time.sleep(3)

	return positive_count

def main(parser):
	options = parser.parse_args()
	edge_ref = options.eref
	fam_ref = options.fref
	fam_id = options.fam
	dir = options.dir	


	if(dir[-1] != "/"):
		dir += "/"

	# create directory for selected pfam
	fam_dir = dir + fam_id + "/"
	os.system("rm -rf %s" %(fam_dir))
	os.system("mkdir -p %s" %(fam_dir))

	# build references
	edge_info = build_edge_guideline(edge_ref)

	positive_count = retrieve_positive(fam_ref, fam_id, edge_info, fam_dir)


	# combine data
	os.system("echo %d > %stmp.txt" %(positive_count, fam_dir))
	os.system("cat %s*_*.txt >> %stmp.txt" %(fam_dir, fam_dir))

	# convert aa to number
	aa_to_number( "%s/tmp.txt" %(fam_dir), "%s/positive_graphs.txt" %(fam_dir), [2])

	echo("There are %d graphs for %s" %(positive_count, fam_id))

	echo("Done")



if __name__ == "__main__":

        parser = argparse.ArgumentParser(prog='03_select_positive_sets.py')
        parser.add_argument("-e", "--edge_reference", dest = "eref", type=str, help="edge reference file", required = True)
       	parser.add_argument("-r", "--family_reference", dest = "fref", type=str, help="pfam or scop association file", required = True)
	parser.add_argument("-f", "--family_id", dest = "fam", type=str, help="family id", required = True)
	parser.add_argument("-d", "--directory", dest = "dir", type=str, help="directory for output", required = True)
	main(parser)
