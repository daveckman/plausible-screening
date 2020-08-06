#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 1

#$ -l h_vmem=16g

#$ -l h_rt=24:00:00

##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; run_PO_tandem_iid_timings($R, $discrep_index, $mode); catch me; display(me); end; quit"

## call with: qsub -v R=1000,discrep_index=1,mode=1 run_script7.sh