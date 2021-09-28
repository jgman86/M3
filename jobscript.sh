#!/bin/bash

#SBATCH -p parallel
#SBATCH -C broadwell
#SBATCH -n 1                    
#SBATCH -c 40
#SBATCH -t 05:30:00              # Run time (hh:mm:ss)
#SBATCH --mem=498000

##SBATCH -A m2_jgu-sim3
##SBATCH -A m2_zdvhpc

module purge # ensures vanilla environment
module load lang/R # will load most current version of R

# for testing
srun Simulation_CS_LKJ.R -N $1 -K $2 --nRetrievals $3 --reps2con $4
