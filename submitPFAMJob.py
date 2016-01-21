import os, time, random, sys, traceback
import subprocess
from subprocess import Popen, PIPE

INPUT = "/u/home/c/chelseaj/project/ProteinGraph/selected_pfam.txt"
OUTDIR = "/u/scratch/c/chelseaj/ProteinGraph/"
SCRIPTDIR = "/u/scratch/c/chelseaj/ProteinGraph/scripts"

TEMPLATE_SERIAL = """
#####################################
#$ -S /bin/bash
#$ -cwd
#$ -N {name}
#$ -e {errfile}
#$ -o {logfile}
#$ -pe make {slots}
#####################################
. /u/local/Modules/default/init/modules.sh
module load python/2.7

cd /u/home/c/chelseaj/project/ProteinGraph
echo $PWD
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""

# create directory for scripts and output
#os.system("mkdir -p %s %s" %(OUTDIR, SCRIPTDIR))

fh = open(INPUT, 'r')
for line in fh:
	(pfam, count) = line.rstrip().split()

	#make the bash file to submit
	scriptfile = SCRIPTDIR + "/" + pfam + ".qsub"
	logfile = SCRIPTDIR + "/" + pfam + ".log"
	errfile = SCRIPTDIR + "/" + pfam + ".err"

	
	script = "python 03_select_proteins.py " + \
		" -e edge_info_v1.txt " + \
		" -r pdb_pfam_mapping.txt " + \
		" -f " + pfam + "_e1" +\
		" -t pfam " + \
		" -d " + OUTDIR + "\n"

        script += "python 03_select_proteins.py " + \
                " -e edge_info_v2.txt " + \
                " -r pdb_pfam_mapping.txt " + \
                " -f " + pfam + "_e2" +\
                " -t pfam " + \
                " -d " + OUTDIR

	scriptFILEHandler = open(scriptfile, 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name=pfam, logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('qsub -cwd -V -l h_data=1G,h_rt=2:00:00 ' + scriptfile,  shell=True)

fh.close()



