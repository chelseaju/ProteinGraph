"""
	Examine the number of proteins and scop domain 
"""

import os, re, argparse
from echo import *

PDB_DIC = {}
SCOP_DIC = {}


def main(parser):
	options = parser.parse_args()
	infile = options.infile
	outdir = options.outdir

	if(outdir[-1] != "/"):
		outdir = outdir + "/"

	# create new directory
	os.system("mkdir -p %s" %(outdir))

	# output file
	pdb_summary = outdir + "pdb_count.txt"
	scop_summary = outdir + "scop_count.txt"

	# retrieve data
	echo("Retrieving data from %s" %(infile))
	fh = open(infile, 'r')
	for line in fh:
		if(line[0] != '#'):
			data = line.split("\t")
			pdb_id = data[1]
			chain = data[2][0]
	
			annotation = dict(item.split("=") for item in data[5].split(",")) 
			scop = annotation["fa"]
			
			pdb_name = pdb_id + "_" + chain

			if(PDB_DIC.has_key(pdb_name)):
				PDB_DIC[pdb_name].append(scop)
			else:
				PDB_DIC[pdb_name] = [scop]
		
			if(SCOP_DIC.has_key(scop)):
				SCOP_DIC[scop].append(pdb_name)
			else:
				SCOP_DIC[scop] = [pdb_name]


	fh.close()

	# write to data
	pdb_fh = open(pdb_summary, 'w')
	for pdb_key in PDB_DIC.keys():
		pdb_value = PDB_DIC[pdb_key]
		pdb_fh.write("%s\t%d\n" %(pdb_key, len(set(pdb_value))))
	pdb_fh.close()
	echo("Writing data to %s" %(pdb_summary))
	
	scop_fh = open(scop_summary, 'w')
	for scop_key in SCOP_DIC.keys():
		scop_value = SCOP_DIC[scop_key]
		scop_fh.write("%s\t%d\n" %(scop_key, len(set(scop_value))))
	scop_fh.close()
	echo("Writing data to %s" %(scop_summary))
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="01_scop_profile.py")
	parser.add_argument("-i", "--input", dest="infile", type=str, help="ex: pdb_scop_mapping.txt", required = True)
	parser.add_argument("-o", "--outdir", dest="outdir", type=str, help="ex: scop_summary", required = True)
	main(parser)
