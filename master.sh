#!/bin/bash

declare -ar N=(1 2 3 4 5)
declare -ar K=(1 2 4 8 20)
nRetrievals=100
nCons=2

# Kann Weg
#declare -ar nCons=(2 4)

# now submit jobs for every permutation
# we keep track:
njobs=0

# guess the account
account=$(sacctmgr -n -s list user $USER format=account%30| grep -v none | head -n1 | tr -d " ")
#account=m2_zdvhpc

for n in "${N[@]}"; do
    for k in "${K[@]}"; do
              for reps2con in $(seq 1 100); do
                jobname="stansim.${n}.${k}.${nr}.${reps2con}"
                slurmout="${jobname}.%j.out"
                if [ ! -e ${slurout} ]; then
                    sbatch -A "$account" -J "$jobname" -o "$slurmout" jobscript.sh "$n" "$k" "$nRetrievals" "$nCons"
                    njobs=$((njobs + 1))
                    if [ $njobs -gt 10000 ]; then
                       exit
                    fi
                fi
            done
        done
    done
