# ProteinGraph

###Purpose:
Convert the protein 3D structure into graph
  
###Description:
Protein graphs are constructed based on the 3D structures. Each vertex denotes an amino acid, and each edge is labeled with the discredited distance between two amino acids. The location of an amino acid is represented by the 3D coordinates of its alpha carbon. Alpha carbon refers to the first carbon atom that attaches to a functional group.

###Data Input:
1. protein and family mapping information: ex pdb_pfam_mapping.txt
2. edge information: define the edge label for different ranges of distance, ex edge_info.txt

###Execution:
* python 01_pfam_profile.py -i pdb_pfam_mapping.txt -o pfam_data_summary
* R --no-save --slave < 02_pfam_distribution.R --args pfam_data_summary/pfam_count.txt pfam_data_summary
* R --no-save --slave < 02_pdb_distribution.R --args pfam_data_summary/pdb_count.txt pfam_data_summary
* python 03_select_protein.py -e edge_info.txt -r pdb_pfam_mapping.txt -f PF01288 -d pfam_graphs

###Details
01: count the number of protein for each family, and number of family for each protein
02: plot the histogram from 01
03: generate the positive and negative graphs for a given family
