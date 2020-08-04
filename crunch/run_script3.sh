#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 4

#$ -l h_vmem=16g

#$ -l h_rt=24:00:00


##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; run_PO_ctsnews_crn($N,$K,$M); catch me; display(me); end; quit"

## call with: qsub -v N=400,K=5,M=3000 run_script3.sh