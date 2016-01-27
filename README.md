# ProteinGraph

###Purpose:
Convert the protein 3D structure into graph
  
###Description:
Protein graphs are constructed based on the 3D structures. Each vertex denotes an amino acid, and each edge is labeled with the discredited distance between two amino acids. The location of an amino acid is represented by the 3D coordinates of its alpha carbon. Alpha carbon refers to the first carbon atom that attaches to a functional group.

###Data Input:
1. protein and family mapping information: ex pdb_pfam_mapping.txt
2. edge information: define the edge label for different ranges of distance, ex edge_info.txt

###Data Input v2:
1. a list of high quality protein filtered by PISCES (filtered_pdb_entries.txt)
2. protein and family mapping information (pdb_pfam_mapping.txt, pdb_scop_mapping.txt)
3. edge information: define the edge label for different ranges of distance (edge_info.txt)

###Scripts:
* 01_pfam_profile.py: examine the number of proteins and families in pfam
* 01_scop_profile.py: examine the number of proteins and families in scop
* 01_assign_family.py: associate high quality protein from PDB with pfam / scop
* 02_pfam_distribution.R: plot the frequencies of proteins and families in pfam
* 02_pfam_distribution.R: plot the frequencies of proteins and families in scop
* 02_family_profile.py: examine the number of proteins in each family (similar to 01_pfam_profile.py, 01_scop_profile.py)
* 02_family_distribution.R: plot the frequencies of proteins of each family
* 03_select_protein.py: select the positive set and negative set (obsolete)
* 03_select_positive_sets.py: generate protein graphs for a given family
* 03_select_negative_sets.py: randomly select from high quality proteins
* 04_extract_embedding.py: extract the amino acid sequence of a given pfam id
* 05_proteins_to_graph.py: convert the proteins in a given family to graph format

###Execution:
* python 01_pfam_profile.py -i pdb_pfam_mapping.txt -o pfam_data_summary
* R --no-save --slave < 02_pfam_distribution.R --args pfam_data_summary/pfam_count.txt pfam_data_summary
* R --no-save --slave < 02_pdb_distribution.R --args pfam_data_summary/pdb_count.txt pfam_data_summary
* python 05_proteins_to_graph.py -e edge_info.txt -r pdb_pfam_mapping.txt -f PF01288 -t pfam -d pfam_dir  
* (obsolete) python 03_select_protein.py -e edge_info.txt -r pdb_pfam_mapping.txt -f PF01288 -d pfam_graphs

###Execution v2:
* python 01_assign_family.py -i filtered_pdb_entries.txt -d pdb_pfam_mapping.txt -t pfam -o filtered_pfam.txt
* python 01_assign_family.py -i filtered_pdb_entries.txt -d pdb_scop_mapping.txt -t scop -o filtered_scop.txt
* python 02_family_profile.py -i filtered_pfam.txt -o filtered_pfam_profile.txt
* python 02_family_profile.py -i filtered_scop.txt -o filtered_scop_profile.txt
* R --no-save --slave < 02_family_distribution.R --args filtered_pfam_profile.txt filtered_pfam_hist.png
* R --no-save --slave < 02_family_distribution.R --args filtered_scop_profile.txt filtered_scop_hist.png
* python 03_select_positive_sets.py -e edge_info_v1.txt -r pdb_pfam_mapping.txt -f PF00001 -d positive_e1
* python 03_select_negative_sets.py -e edge_info_v1.txt -p filtered_pdb_entries.txt -n 20 -d negative_e1

###Sources:
* pfam - ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/
* SCOP - http://scop.berkeley.edu/downloads/parse/dir.cla.scope.2.05-stable.txt
* PISCES - http://dunbrack.fccc.edu/PISCES.php
