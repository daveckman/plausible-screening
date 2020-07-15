#!/bin/bash

#$ -cwd

#$ -j y

#$ -S /bin/bash

#$ -pe smp 4

#$ -l h_vmem=400g


##export OMP_NUM_THREADS=1
/share/apps/matlab_current/bin/matlab -r "try; run_PO_sS; catch me; display(me); end; quit"
##/share/apps/matlab_current/bin/matlab -r "try; run_PO_sS($it, $N, $p); catch me; display(me); end; quit"