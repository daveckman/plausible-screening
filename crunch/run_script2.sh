#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 8

#$ -l h_vmem=16g

#$ -l h_rt=24:00:00

##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; run_PO_tandem_iid($m, $discrep_index, $mode); catch me; display(me); end; quit"

## call with: qsub -v m=1,discrep_index=3,mode=1 run_script2.sh