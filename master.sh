#!/bin/bash

#SBATCH -p parallel
#SBATCH -C broadwell
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -t 05:30:00              # Run time (hh:mm:ss)


#SBATCH -A m2_jgu-sim3

declare -ar N=(1 2)
declare -ar K=(8 20)
declare -ar nCon=(2 4)
declare -ar nRetrievals=(100)


# now submit jobs for every permutation
# we keep track:
njobs=0

# guess the account
account=$(sacctmgr -n -s list user $USER format=account%30| grep -v none | head -n1 | tr -d " ")
#account=m2_zdvhpc

for n in "${N[@]}"; do
    for k in "${K[@]}"; do
      for j in "${nCons[@]}"; do
              for reps2con in $(seq 1 10); do
                jobname="stansim.${n}.${k}.${nr}.${nCon}.${reps2con}"
                slurmout="${jobname}.%j.out"
                if [ ! -e ${slurout} ]; then
                    sbatch -A "$account" -J "$jobname" -o "$slurmout" jobscript.sh "$n" "$k" "$nRetrievals" "$nCon"
                    njobs=$((njobs + 1))
                    if [ $njobs -gt 10000 ]; then
                       exit
                    fi
                fi
            done
        done
    done
done
    
