#!/bin/bash
#PBS -N testrun
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2gb
#PBS -l walltime=00:02:00
#PBS -j oe
#PBS -q gpu
#PBS -M "emailaddress@somewhere.xx"
#PBS -m ae

# some debug output
echo $PBS_O_WORKDIR
hostname

cd $PBS_O_WORKDIR

module load atistream/2.3

./hmc.v0.1 input_heatbath_benchmark 0
