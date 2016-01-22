# ProteinGraph

###Purpose:
Convert the protein 3D structure into graph
  
###Description:
Protein graphs are constructed based on the 3D structures. Each vertex denotes an amino acid, and each edge is labeled with the discredited distance between two amino acids. The location of an amino acid is represented by the 3D coordinates of its alpha carbon. Alpha carbon refers to the first carbon atom that attaches to a functional group.

###Data Input:
1. protein and family mapping information: ex pdb_pfam_mapping.txt
2. edge information: define the edge label for different ranges of distance, ex edge_info.txt

###Scripts:
* 01_pfam_profile.py: examine the number of proteins and families in pfam
* 01_scop_profile.py: examine the number of proteins and families in scop
* 02_pfam_distribution.R: plot the frequencies of proteins and families in pfam
* 02_pfam_distribution.R: plot the frequencies of proteins and families in scop
* 03_select_protein.py: select the positive set and negative set (obsolete)
* 04_extract_embedding.py: extract the amino acid sequence of a given pfam id
* 05_proteins_to_graph.py: convert the proteins in a given family to graph format

###Execution:
* python 01_pfam_profile.py -i pdb_pfam_mapping.txt -o pfam_data_summary
* R --no-save --slave < 02_pfam_distribution.R --args pfam_data_summary/pfam_count.txt pfam_data_summary
* R --no-save --slave < 02_pdb_distribution.R --args pfam_data_summary/pdb_count.txt pfam_data_summary
* python 05_proteins_to_graph.py -e edge_info.txt -r pdb_pfam_mapping.txt -f PF01288 -t pfam -d pfam_dir  
* (obsolete) python 03_select_protein.py -e edge_info.txt -r pdb_pfam_mapping.txt -f PF01288 -d pfam_graphs

###Sources:
* pfam - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/
* SCOP - http://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.05-stable.txt
