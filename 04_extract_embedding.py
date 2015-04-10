"""
	Extract the embedding amino acid sequence
"""
import os, argparse, subprocess
import Bio.PDB.Polypeptide
import echo

def retrieve_fasta(infile, outfile, pfam):

	echo.echo("Reading data from %s" %(infile))
	out_fh = open(outfile, 'w')
	fh = open(infile, 'r')
	for line in fh:
		line = line.rstrip()
		if(line[0] == 'G'):
			details = line.split()
			graph_id = details[0][2:]
			nodes = details[1][2:].split(',')
			nodes = [int(x) for x in nodes]
			nodes = sorted(nodes)

			p1 = subprocess.Popen(['grep', 't # %s ' %(graph_id), "pfam_graphs/%s/positive/positive_graphs.txt" %(pfam)],  stdout= subprocess.PIPE)
			protein_id = p1.communicate()[0].split()[-1]
			seq = ""
			index = 1
			for n in nodes:
				p2 = subprocess.Popen(['grep', '^v %s ' %(n), "pfam_graphs/%s/positive/%s.txt" %(pfam, protein_id)], stdout = subprocess.PIPE)
				aa = p2.communicate()[0].split()[-1]
				
				try:
					seq_a = Bio.PDB.Polypeptide.three_to_one(aa)
				except:
					seq_a = "?"


				while(index < int(n)):
					seq += "*"
					index += 1
				
				seq += seq_a
				index += 1
		
			out_fh.write(">%s\n" %(protein_id))
			out_fh.write("%s\n" %(seq))

	fh.close()
	out_fh.close()

	echo.echo("Writing sequence into %s" %(outfile))


def main(parser):
	options = parser.parse_args()
	dir = options.dir
	pfam = options.pfam

	if(dir[-1] != "/"):
		dir = dir + "/"
	
	dir = dir + pfam
	
	for p_file in os.listdir(dir):
		if(p_file.startswith("pattern")):
			infile = dir + "/" + p_file
			outfile = dir + "/" + "seq_" + p_file + ".fa"
	
			retrieve_fasta(infile, outfile, pfam)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='04_extract_embeding.py')
	parser.add_argument("-d", "--dir", dest = "dir", type=str, help="input and output directory", required = True)
	parser.add_argument("-f", "--pfam", dest = "pfam", type=str, help="pfam", required = True)
	main(parser)
