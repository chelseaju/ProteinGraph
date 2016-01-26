"""
	For a list of high qualities proteins, link their pfam id and scop id
"""

import os, re, argparse
from echo import *

PDB = {}

"""
	Read in high quality PDB chains
"""
def build_pdb(pdb):
	fh = open(pdb, 'r')
	for line in fh:
		pdb_id = line.split()[0]
		PDB[pdb_id] = ""
	fh.close()

"""
	Associate with pfam
"""
def build_pfam_db(pfam):
	fh = open(pfam, 'r')
	for line in fh:
		data = line.split()
		pdb_id = data[0] + data[1]
		pfam_id = data[4]
		print pfam_id

	fh.close()


def main(parser):
	options = parser.parse_args()
	pdb = options.pdb
	pfam = options.db
	db_type = options.dbtype
	outfile = options.outfile

	build_pdb(pdb)
	if(db_type == "pfam"):
		build_pfam_db(pfam)
		
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(prog="06_assign_family.py")	
	parser.add_argument("-i", "--pdb", dest = "pdb", type=str, help="high quality pdb", required = True)
	parser.add_argument("-d", "--database", dest="db", type=str, help="mapping information", required = True)
	parser.add_argument("-t", "--type", dest="dbtype", type=str, help="type of maaping, ie pfam, scop", required = True)
	parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="output filename", required = True)
	main(parser)
