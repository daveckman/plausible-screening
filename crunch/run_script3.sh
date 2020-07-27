#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 4

#$ -l h_vmem=16g

#$ -l h_rt=24:00:00


##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; run_PO_ctsnews_crn($K,$N); catch me; display(me); end; quit"

## call with: qsub -v K=20,N=400 run_script3.sh