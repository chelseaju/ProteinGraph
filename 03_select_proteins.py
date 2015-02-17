"""
	Select positive and negative sets of protein for a given family
"""


import argparse, re, os, time, random
from pdb import *
from graph import *
from echo import *


POSITIVE = []
NEGATIVE = {}

# retrieve the list of potential positive and negative proteins
def build_pfam_reference(pfam_ref, pfam_id):

	fh = open(pfam_ref, 'r')
        for line in fh:
                data = line.split("\t")
                pdb_id = data[0]
		chain = data[1]
                pfam = data[4]

		if(not pdb_id == "PDB_ID"):
			pdb_name = pdb_id + "_" + chain
			pfam_match = re.match(r'(.*)\.(\d+)', pfam)
                        if(pfam_match):
                                pfam = pfam_match.group(1)
				if(pfam == pfam_id):
					POSITIVE.append(pdb_name)
				else:
					if(NEGATIVE.has_key(pdb_name)):
						NEGATIVE[pdb_name].append(pfam)
					else:
						NEGATIVE[pdb_name] = [pfam]

        fh.close()
	
	echo("Building Pfam Reference")


# positive graphs
def retrieve_positive(pfam_id, edge_info, dir):
	positive_count = 0
	for pdb_chain in sorted(set(POSITIVE)):
		(pdb, chain) = pdb_chain.split('_')
		
		# distance file name
		pdb_dst_file = dir + "/" + pdb_chain + ".dist"

		# graph file name
		pdb_graph_file = dir + "/" + pdb_chain + ".txt"

		# parse pdb data
		pdb_info = parse_pdb_by_id(pdb, chain)

		if(pdb_info):
			# comput distance
			pw_dist = pairwise_distance(pdb_info, pdb_dst_file)
			
			# convert structure to graph
			title = pfam_id+ " " + pdb_chain
			pdb_to_graph(pw_dist, pdb_graph_file, edge_info, positive_count, title)
		
			positive_count = positive_count + 1

		if(positive_count % 100 == 0):
			time.sleep(3)

	return positive_count

# negative graphs
def retrieve_negative(pfam_id, edge_info, dir, count):

	# negative candidates
	candidates = NEGATIVE.keys()
	selected = []

	random.seed()

	count = count * (1.5)
	negative_count = 0
	while(negative_count < count):
		index = random.randint(0, len(candidates)-1)
		negative_protein = candidates[index]
	
		if(not negative_protein in selected and not negative_protein in POSITIVE):
			selected.append(negative_protein)
			negative_count = negative_count + 1
		

	negative_count = 0
	for pdb_chain in sorted(selected):
		(pdb, chain) = pdb_chain.split('_')

		# distance file name
		pdb_dst_file = dir + "/" + pdb_chain + ".dist"
		
		# graph file name
		pdb_graph_file = dir + "/" + pdb_chain + ".txt"

		# parse pdb data
		pdb_info = parse_pdb_by_id(pdb, chain)
	
                if(pdb_info):
                        # comput distance
                        pw_dist = pairwise_distance(pdb_info, pdb_dst_file)

                        # convert structure to graph
                        title = pdb_chain + " " + ",".join(set(NEGATIVE[pdb_chain]))
                        pdb_to_graph(pw_dist, pdb_graph_file, edge_info, negative_count, title)

                        negative_count = negative_count + 1

                if(negative_count % 100 == 0):
                        time.sleep(3)


	return negative_count
		

def main(parser):
	options = parser.parse_args()
	edge_ref = options.eref
	pfam_ref = options.pref
	pfam_id = options.pfam
	dir = options.dir
	
	if(dir[-1] != "/"):
		dir += "/"

	# create directory for selected pfam
	pfam_dir = dir + pfam_id + "/"
	pfam_positive_dir = pfam_dir + "positive"
	pfam_negative_dir = pfam_dir + "negative"
	os.system("rm -rf %s" %(pfam_negative_dir))
	os.system("rm -rf %s" %(pfam_positive_dir))
	os.system("mkdir -p %s %s %s" %(pfam_dir, pfam_positive_dir, pfam_negative_dir))

	# build references
	edge_info = build_edge_guideline(edge_ref)
	build_pfam_reference(pfam_ref, pfam_id)


	positive_count = retrieve_positive(pfam_id, edge_info, pfam_positive_dir)
	negative_count = retrieve_negative(pfam_id, edge_info, pfam_negative_dir, positive_count)


	# merge negative and positive graph
	os.system("echo %d > %s/tmp.txt" %(positive_count, pfam_positive_dir))
	os.system("echo %d > %s/tmp.txt" %(negative_count, pfam_negative_dir))
	os.system("cat %s/*_*.txt >> %s/tmp.txt" %(pfam_positive_dir, pfam_positive_dir))
	os.system("cat %s/*_*.txt >> %s/tmp.txt" %(pfam_negative_dir, pfam_negative_dir))

	aa_to_number( "%s/tmp.txt" %(pfam_positive_dir), "%s/positive_graphs.txt" %(pfam_positive_dir), [2])
        aa_to_number( "%s/tmp.txt" %(pfam_negative_dir), "%s/negative_graphs.txt" %(pfam_negative_dir), [2])

	echo("There are %d positive and %d negative graphs for %s" %(positive_count, negative_count, pfam_id))
	echo("Writing positive graph to %s/positive_graphs.txt" %(pfam_positive_dir))
	echo("Writing negative graph to %s/negative_graphs.txt" %(pfam_negative_dir))

	echo("Done")
if __name__ == "__main__":

        parser = argparse.ArgumentParser(prog='03_select_protein.py')
        parser.add_argument("-e", "--edge_reference", dest = "eref", type=str, help="edge reference file", required = True)
       	parser.add_argument("-r", "--pfam_reference", dest = "pref", type=str, help="pfam association file", required = True)
	parser.add_argument("-f", "--pfam_id", dest = "pfam", type=str, help="pfam id", required = True)
	parser.add_argument("-d", "--directory", dest = "dir", type=str, help="directory for output", required = True)
	main(parser)
