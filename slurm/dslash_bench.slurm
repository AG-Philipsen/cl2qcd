#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=OpTiMaL_dslash_bench.log
#SBATCH --job-name=OpTiMaL_dslash_bench
#SBATCH --mem-per-cpu=400
#SBATCH --mail-type=ALL
#SBATCH --partition=parallel
#SBATCH --gres=gpu
#SBATCH --time=0-02:00:00
#SBATCH --dependency=singleton
#

export DISPLAY=:0

# report where we are
srun hostname

# check gpu
srun aticonfig --odgc --odgt --adapter=all

# blablabla
srun aticonfig --od-enable
# modify gpu clock to factory defaults
srun aticonfig --odsc 850,1200

# check gpu again
srun aticonfig --odgc --odgt --adapter=all

echo "---------------------------"

srun ../../benchmarks/dslash_bench.py

echo "---------------------------"

# Reset clocks to default
srun aticonfig --odsc 750,1100

# done messing with clocks
srun aticonfig --od-disable

# check gpu again
srun aticonfig --odgc --odgt --adapter=all
