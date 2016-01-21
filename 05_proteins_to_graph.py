"""
	Retrieve proteins from a given family
	Convert these proteins into graph format 
"""


import argparse, re, os, time, random
from pdb import *
from graph import *
from echo import *

POSITIVE = []
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
        fh.close()

        echo("Building Pfam Reference")

# positive graphs
def retrieve_positive(pfam_id, edge_info, dir):
        positive_count = 0
        positive_candidates = set(POSITIVE)
        select_candidates = random.sample(positive_candidates, 11)

        for pdb_chain in sorted(select_candidates):
                (pdb, chain) = pdb_chain.split('_')

                # distance file name
                pdb_dst_file = dir + pdb_chain + ".dist"

                # graph file name
                pdb_graph_file = dir + pdb_chain + ".txt"

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

# positive graphs
def retrieve_graph(pfam_id, edge_info, dir):
        positive_count = 0
        positive_candidates = set(POSITIVE)
     #   select_candidates = random.sample(positive_candidates, 11)

        for pdb_chain in sorted(positive_candidates):
                (pdb, chain) = pdb_chain.split('_')

                # distance file name
                pdb_dst_file = dir + pdb_chain + ".dist"

                # graph file name
                pdb_graph_file = dir + pdb_chain + ".txt"

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
	os.system("mkdir -p %s " %(fam_dir))

	# build references
        edge_info = build_edge_guideline(edge_ref)

        if(ftype == "pfam"):
                build_pfam_reference(fam_ref, fam_id)
        elif(ftype == "scop"):
                build_scop_reference(fam_ref, fam_id)

	num_protein = retrieve_graph(fam_id, edge_info, fam_dir)
	echo("Retrieving %d proteins for %s" %(num_protein, fam_id)) 

if __name__ == "__main__":

        parser = argparse.ArgumentParser(prog='05_proteins_to_graph.py')
        parser.add_argument("-e", "--edge_reference", dest = "eref", type=str, help="edge reference file", required = True)
        parser.add_argument("-r", "--family_reference", dest = "fref", type=str, help="pfam or scop association file", required = True)
        parser.add_argument("-f", "--family_id", dest = "fam", type=str, help="family id", required = True)
        parser.add_argument("-t", "--family_type", dest = "ftype", type=str, help="pfam or scop", required = True)
        parser.add_argument("-d", "--directory", dest = "dir", type=str, help="directory for output", required = True)
        main(parser)
