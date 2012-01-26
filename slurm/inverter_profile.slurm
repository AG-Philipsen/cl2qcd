#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=OpTiMaL_profile.log
#SBATCH --job-name=OpTiMaL_profile
#SBATCH --mem-per-cpu=400
#SBATCH --mail-type=ALL
#SBATCH --partition=parallel
#SBATCH --gres=gpu
#SBATCH --time=0-00:15:00
#SBATCH --dependency=singleton
#

export DISPLAY=:0

# report where we are
srun hostname

# check gpu
srun aticonfig --odgc --odgt --adapter=all

# run inverter
srun ./inverter_benchmarks ../inverter_bench_input

# check gpu again
srun aticonfig --odgc --odgt --adapter=all
