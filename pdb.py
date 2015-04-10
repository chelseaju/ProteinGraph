"""
	Retrieve protein from PDB
"""

from Bio.PDB import *
from echo import *

PDBL = PDBList()
PDBP = PDBParser()

"""
	Convert 3 letter amino acid code to ascii integer
"""

def encode(aa):
	return ord(aa[0])*10000 + ord(aa[1]) * 100 + ord(aa[2])


"""
	Convert ascii integer to 3 letter amino acid code
"""
def decode(ascii):
	return chr(ascii / 10000) + chr((ascii % 10000) / 100) + chr(((ascii % 10000) % 100))


"""
	Download PDB and parse the alpha carbon information
"""
def parse_pdb_by_id(pdb_id, chain_label):
	
	echo("Parsing PDB %s" %(pdb_id))
	try:
		pdb_file = PDBL.retrieve_pdb_file(pdb_id)
	except:
		return None

	chain = PDBP.get_structure(pdb_id, pdb_file)[0][chain_label]

	seq = []

	residue_count = 1
	for residue in chain:
		if(is_aa(residue) and residue.has_id('CA')):
			residue_id = residue.get_resname() + "_" + str(residue_count)
			residue_count = residue_count + 1

			seq.append((residue_id, residue['CA']))
	return seq

"""
	Comput pairwise distance between two amino acid based on alpha carbon
"""
def pairwise_distance(pdb_info, pdb_dst_file):

	distance = []
	fh = open(pdb_dst_file, 'w')
	for i in xrange(len(pdb_info)):
		for j in xrange(i+1, len(pdb_info)):
			distance.append((pdb_info[i][0], pdb_info[j][0], pdb_info[i][1] - pdb_info[j][1]))
			fh.write("%s %s %f\n" %(pdb_info[i][0], pdb_info[j][0], pdb_info[i][1] - pdb_info[j][1]))
	fh.close()

	echo("Writing distance to %s" %(pdb_dst_file))

	return distance

"""
	Convert AA symbol to number
"""
def aa_to_number(infile, outfile, columns):
	infh = open(infile, 'r')
	outfh = open(outfile, 'w')

	for line in infh:
		if(line[0] == 'v'):
			data = line.rstrip().split()

			for c in columns:
				data[c] = str(encode(data[c]))
		
			outfh.write("%s\n" %(" ".join(data)))
		else:
			outfh.write(line)
	infh.close()
	outfh.close()
		

if __name__ == "__main__":

#	pdb_info = parse_pdb_by_id('3CKZ', 'A')
	pdb_info = parse_pdb_by_id('1L2Y', 'A')
#	pdb_info = parse_pdb_by_id('4BPP', 'A')
	pdb_dist = pairwise_distance(pdb_info, "test.txt")

#	ascii = encode('MET')
#	print ascii
#	print decode(ascii)


#	aa_to_number("pfam_graphs/PF03645/positive/tmp.txt", "pfam_graphs/PF03645/positive/positive_graphs.txt", [2])


