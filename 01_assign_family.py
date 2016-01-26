"""
	For a list of high qualities proteins, link their pfam id and scop id
"""

import os, re, argparse, datetime

PDB = {}

def echo(message):
        print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(message))


"""
	Read in high quality PDB chains
"""
def build_pdb(pdb):
	fh = open(pdb, 'r')
	for line in fh:
		pdb_id = line.split()[0]
		PDB[pdb_id] = []
	fh.close()

"""
	Associate with pfam
"""
def build_pfam_db(pfam):
	fh = open(pfam, 'r')
	for line in fh:
		data = line.split()
		pdb_id = data[0] + data[1]
		pfam_id = data[4][0:7]
		
		if(PDB.has_key(pdb_id)):
			PDB[pdb_id].append(pfam_id)

	fh.close()

"""
	Associate with scop
"""
def build_scop_db(scop):
	fh = open(scop, 'r')
	for line in fh:
		if(line[0] != '#'):
			data = line.split()
			pdb_id = data[0][1:6].upper()
			scop_annotation = dict(item.split("=") for item in data[5].split(","))
			scop_id = scop_annotation["fa"]
			
			if(PDB.has_key(pdb_id)):
				PDB[pdb_id].append(scop_id)
	fh.close()


def export_data(outfile):
	fh = open(outfile, 'w')
	for k in PDB.keys():
		for fam in set(PDB[k]):
			fh.write("%s\t%s\n" %(k, fam))
	fh.close()

def main(parser):
	options = parser.parse_args()
	pdb = options.pdb
	db = options.db
	db_type = options.dbtype
	outfile = options.outfile

	echo("Reading Filtered PDB")
	build_pdb(pdb)
	if(db_type == "pfam"):
		echo("Pfam Association")
		build_pfam_db(db)
		
	elif(db_type =="scop"):
		echo("Scop Association")
		build_scop_db(db)

	echo("Writing to file")
	export_data(outfile)	

	echo("Done")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(prog="01_assign_family.py")	
	parser.add_argument("-i", "--pdb", dest = "pdb", type=str, help="high quality pdb", required = True)
	parser.add_argument("-d", "--database", dest="db", type=str, help="mapping information", required = True)
	parser.add_argument("-t", "--type", dest="dbtype", type=str, help="type of maaping, ie pfam, scop", required = True)
	parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="output filename", required = True)
	main(parser)
