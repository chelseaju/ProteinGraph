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


def build_scop_reference(scop_ref, fam_id):
	fh = open(scop_ref, 'r')
	for line in fh:
		if(line[0] != "#"):
			data = line.split("\t")
                        pdb_id = data[1]
                        chain = data[2][0]

                        annotation = dict(item.split("=") for item in data[5].split(","))
                        scop = annotation["fa"]

                        pdb_name = pdb_id + "_" + chain

			if(fam_id == scop):
				POSITIVE.append(pdb_name)
			else:
				if(NEGATIVE.has_key(pdb_name)):
					NEGATIVE[pdb_name].append(scop)
				else:
					NEGATIVE[pdb_name] = [scop]
	fh.close()
	echo("Building Scop Reference")



# positive graphs
def retrieve_positive(pfam_id, edge_info, dir):
	positive_count = 0
	positive_candidates = set(POSITIVE)
	select_candidates = random.sample(positive_candidates, 11)

	for pdb_chain in sorted(select_candidates):
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

	count = count * 1.5
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
	fam_ref = options.fref
	fam_id = options.fam
	ftype = options.ftype
	dir = options.dir	


	if(dir[-1] != "/"):
		dir += "/"

	# create directory for selected pfam
	fam_dir = dir + fam_id + "/"
	fam_positive_dir = fam_dir + "positive"
	fam_negative_dir = fam_dir + "negative"
	os.system("rm -rf %s" %(fam_negative_dir))
	os.system("rm -rf %s" %(fam_positive_dir))
	os.system("mkdir -p %s %s %s" %(fam_dir, fam_positive_dir, fam_negative_dir))

	# build references
	edge_info = build_edge_guideline(edge_ref)

	if(ftype == "pfam"):
		build_pfam_reference(fam_ref, fam_id)
	elif(ftype == "scop"):
		build_scop_reference(fam_ref, fam_id)

	positive_count = retrieve_positive(fam_id, edge_info, fam_positive_dir)
	negative_count = retrieve_negative(fam_id, edge_info, fam_negative_dir, positive_count)


	# merge negative and positive graph
	os.system("echo %d > %s/tmp.txt" %(positive_count, fam_positive_dir))
	os.system("echo %d > %s/tmp.txt" %(negative_count, fam_negative_dir))
	os.system("cat %s/*_*.txt >> %s/tmp.txt" %(fam_positive_dir, fam_positive_dir))
	os.system("cat %s/*_*.txt >> %s/tmp.txt" %(fam_negative_dir, fam_negative_dir))

	aa_to_number( "%s/tmp.txt" %(fam_positive_dir), "%s/positive_graphs.txt" %(fam_positive_dir), [2])
        aa_to_number( "%s/tmp.txt" %(fam_negative_dir), "%s/negative_graphs.txt" %(fam_negative_dir), [2])

	echo("There are %d positive and %d negative graphs for %s" %(positive_count, negative_count, fam_id))
	echo("Writing positive graph to %s/positive_graphs.txt" %(fam_positive_dir))
	echo("Writing negative graph to %s/negative_graphs.txt" %(fam_negative_dir))

	echo("Done")
if __name__ == "__main__":

        parser = argparse.ArgumentParser(prog='03_select_protein.py')
        parser.add_argument("-e", "--edge_reference", dest = "eref", type=str, help="edge reference file", required = True)
       	parser.add_argument("-r", "--family_reference", dest = "fref", type=str, help="pfam or scop association file", required = True)
	parser.add_argument("-f", "--family_id", dest = "fam", type=str, help="family id", required = True)
	parser.add_argument("-t", "--family_type", dest = "ftype", type=str, help="pfam or scop", required = True)
	parser.add_argument("-d", "--directory", dest = "dir", type=str, help="directory for output", required = True)
	main(parser)
