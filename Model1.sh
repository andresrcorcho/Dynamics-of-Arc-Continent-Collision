#!/bin/bash
#PBS -P m18
#PBS -q express
#PBS -l walltime=01:59:00
#PBS -l mem=128GB
#PBS -l jobfs=1GB
#PBS -l ncpus=32
## For licensed software, you have to specify it to get the job running. For unlicensed software, you should also specify$
#PBS -l software=my_program
## The job will be executed from current working directory instead of home.
#PBS -l wd

module purge
module load pbs dot gcc/5.2.0 hdf5/1.10.2p petsc/3.8.4 swig/3.0.12 python/2.7.11 openmpi/3.1.2
source /home/563/ac9890/python-2.7.11-venv/bin/activate

MODELNAME="SubductionFirst"
OUTPUTPATH=`pwd`
SCRIPT="MyFirstModel.py"

mpirun -np 32 --mca mpi_warn_on_fork 0 --mca opal_abort_print_stack 1 --mca mpi_param_check 1 \
--mca mpi_add_procs_cutoff 256 python ./$SCRIPT 1> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.log 2> $OUTPUTPATH/$MODELNAME.$PBS_JOBID.err

