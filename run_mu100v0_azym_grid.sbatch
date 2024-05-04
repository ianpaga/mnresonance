#!/bin/bash

#SBATCH --job-name=mnres3
#SBATCH --ntasks=1
#SBATCH --partition=astro2_long
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=32768
#SBATCH --cpus-per-task=20   

echo
echo "SLURM ENVIRONMENT"
echo "-----------------"
env | grep SLURM_ | sort | while read var; do
    echo " * $var"
done

### specify in the prompt which start_ code to submit  

cd start_$1/

EXE="/groups/astro/ianpaga/mnresonance/other_runs/start_$1/mu100-400by400by100v0_azym_grid.sh"

srun time -p $EXE
