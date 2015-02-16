"""
	Examine the number of proteins and pfam
"""

import os, re, argparse
from echo import *

PDB_DIC = {}
PFAM_DIC = {}


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
	pfam_summary = outdir + "pfam_count.txt"

	# retrieve data
	echo("Retrieving data from %s" %(infile))
	fh = open(infile, 'r')
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
				
				if(PFAM_DIC.has_key(pfam)):
					PFAM_DIC[pfam].append(pdb_name)
				else:
					PFAM_DIC[pfam] = [pdb_name]
				
				if(PDB_DIC.has_key(pdb_name)):
					PDB_DIC[pdb_name].append(pfam)
				else:
					PDB_DIC[pdb_name] = [pfam]

	fh.close()

	# write to data
	pdb_fh = open(pdb_summary, 'w')
	for pdb_key in PDB_DIC.keys():
		pdb_value = PDB_DIC[pdb_key]
		pdb_fh.write("%s\t%d\n" %(pdb_key, len(set(pdb_value))))
	pdb_fh.close()
	echo("Writing data to %s" %(pdb_summary))
	
	pfam_fh = open(pfam_summary, 'w')
	for pfam_key in PFAM_DIC.keys():
		pfam_value = PFAM_DIC[pfam_key]
		pfam_fh.write("%s\t%d\n" %(pfam_key, len(set(pfam_value))))
	pfam_fh.close()
	echo("Writing data to %s" %(pfam_summary))
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="01_pfam_profile.py")
	parser.add_argument("-i", "--input", dest="infile", type=str, help="ex: pdb_pfam_mapping.txt", required = True)
	parser.add_argument("-o", "--outdir", dest="outdir", type=str, help="ex: pfam_summary", required = True)
	main(parser)
