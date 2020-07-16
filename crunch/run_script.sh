#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 4

#$ -l h_vmem=400g


##export OMP_NUM_THREADS=1
##/share/apps/matlab_current/bin/matlab -r "try; run_PO_sS; catch me; display(me); end; quit"
/share/apps/matlab_current/bin/matlab -r "try; run_PO_ctsnews_iid($K,$N); catch me; display(me); end; quit"

## call with: qsub -v K=20,N=400 run_script.sh