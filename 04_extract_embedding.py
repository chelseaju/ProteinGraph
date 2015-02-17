"""
	Extract the embedding amino acid sequence
"""
import argparse, subprocess
import Bio.PDB.Polypeptide
import echo

def main(parser):
	options = parser.parse_args()
	infile = options.input
	outfile = options.output
	pfam = options.pfam

	out_fh = open(outfile, 'w')

	fh = open(infile, 'r')
	for line in fh:
		line = line.rstrip()
		if(line[0] == 'g'):
			details = line.split()
			graph_id = details[1]
			nodes = details[2].split(',')

			p1 = subprocess.Popen(['grep', 't # %s' %(graph_id), "pfam_graphs/%s/positive/positive_graphs.txt" %(pfam)],  stdout= subprocess.PIPE)
			protein_id = p1.communicate()[0].split()[-1]
			seq = ""
			index = 1
			for n in nodes:
				p2 = subprocess.Popen(['grep', '^v %s' %(n), "pfam_graphs/%s/positive/%s.txt" %(pfam, protein_id)], stdout = subprocess.PIPE)
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



if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='04_extract_embedding.py')
	parser.add_argument("-i", "--input", dest = "input", type=str, help="embedding subgraphs", required = True)
	parser.add_argument("-o", "--output", dest = "output", type=str, help="output file", required = True)

	parser.add_argument("-f", "--pfam", dest = "pfam", type=str, help="pfam id", required = True)
	main(parser)
