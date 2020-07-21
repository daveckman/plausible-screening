#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 4

#$ -l h_vmem=400g

#$ -l h_rt=24:00:00

##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; run_PO_tandem_iid($M); catch me; display(me); end; quit"

## call with: qsub -v M=10 run_script2.sh