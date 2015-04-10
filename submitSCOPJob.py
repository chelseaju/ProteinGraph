import os, time, random, sys, traceback
import subprocess
from subprocess import Popen, PIPE


INPUT = "selected_scop.txt"
OUTDIR = "scop_graphs"
SCRIPTDIR = "scop_scripts"

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
os.system("mkdir -p %s %s" %(OUTDIR, SCRIPTDIR))

fh = open(INPUT, 'r')
for line in fh:
	(scop, count) = line.rstrip().split()

	#make the bash file to submit
	scriptfile = SCRIPTDIR + "/" + scop + ".qsub"
	logfile = SCRIPTDIR + "/" + scop + ".log"
	errfile = SCRIPTDIR + "/" + scop + ".err"

	
	script = "python 03_select_proteins.py " + \
		" -e edge_info_v2.txt " + \
		" -r pdb_scop_mapping.txt " + \
		" -f " + scop + \
		" -t scop " + \
		" -d " + OUTDIR

	scriptFILEHandler = open(scriptfile, 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="scop_%s"%(scop), logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('qsub -cwd -V -l h_data=1G,h_rt=2:00:00 ' + scriptfile,  shell=True)

fh.close()



