#!/bin/bash

#SBATCH -p parallel
#SBATCH -C broadwell
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -t 05:30:00              # Run time (hh:mm:ss)
#SBATCH -J M3_CS_EE

#SBATCH -A m2_jgu-sim3

module purge # ensures vanilla environment
module load lang/R # will load most current version of R

# for testing
srun  R --vanilla -f Simulation_CS_EE_LKJ.R 
