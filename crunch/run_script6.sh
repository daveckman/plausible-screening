#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 8

#$ -l h_vmem=16g

#$ -l h_rt=24:00:00

##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; fix_linprog_feasibility_errors; catch me; display(me); end; quit"

## call with: qsub run_script6.sh