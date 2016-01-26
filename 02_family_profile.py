"""
	Examine the number of proteins and pfam
"""

import os, re, argparse, datetime

def echo(message):
        print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(message))


FAM_DIC = {}

def main(parser):
	options = parser.parse_args()
	infile = options.infile
	outfile = options.outfile

	# retrieve data
	echo("Retrieving data from %s" %(infile))
	fh = open(infile, 'r')
	for line in fh:
		(pdb_id, fam_id) = line.rstrip().split("\t")
		
		if(FAM_DIC.has_key(fam_id)):
			FAM_DIC[fam_id] += 1
		else:
			FAM_DIC[fam_id] = 1	
	
	fh.close()

	# write to data
	echo("Writing summary to %s" %(outfile))
	outfh = open(outfile, 'w')
	for fam in FAM_DIC.keys():
		outfh.write("%s\t%d\n" %(fam, FAM_DIC[fam]))
	outfh.close()
	
	echo("Done")	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="02_family_profile.py")
	parser.add_argument("-i", "--input", dest="infile", type=str, help="ex: pfam_filtered.txt", required = True)
	parser.add_argument("-o", "--outfile", dest="outfile", type=str, help="ex: pfam_summary.txt", required = True)
	main(parser)
