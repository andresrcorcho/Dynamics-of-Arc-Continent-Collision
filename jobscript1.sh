#!/bin/bash
#PBS -P m18
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=128GB
#PBS -l jobfs=1GB
#PBS -l ncpus=16
## For licensed software, you have to specify it to get the job running. For unlicensed software, you should also specify$
#PBS -l software=my_program
## The job will be executed from current working directory instead of home.
#PBS -l wd


module load underworld/2.5.0
mpirun -np 16 python localscript1.py > my_output1.out
