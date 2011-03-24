# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil
import loewe_defs

#the walltime setting have to be adjusted!!!
jobscript = """
#!/bin/bash
#PBS -N heatbath_benchmark
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2gb
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -q """

jobscript2 = """
#PBS -M "pinke@th.physik.uni-frankfurt.de"
#PBS -m ae

# some debug output
echo $PBS_O_WORKDIR
hostname

cd $PBS_O_WORKDIR

module load atistream/2.3

./hmc ../../input_heatbath_benchmark """

globaldefs1 = """
#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

#ifndef _INKERNEL_

"""

globaldefs2 = """
#define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

#define SPINORSIZE NSPIN*NC
#define HALFSPINORSIZE NSPIN/2*NC
#define SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE
#define EOPREC_SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE/2

//startconditions:
#define START_FROM_SOURCE 2
#define COLD_START 0
#define HOT_START 1

#endif //_INKERNEL_

//EVEN ODD
#define EVEN 0
#define ODD 1

#define TRUE 1
#define FALSE 0

#define PI 	3.14159265358979

#define su2_entries 4

#ifdef _USEGPU_
#define NUMTHREADS 128
#else
#define NUMTHREADS 1
#endif

#endif
"""

cur_dir = os.getcwd()

#suppose there is a dir for every option-setting (8 dir) with a suiting build dir in it
#this can be changed via cmake -DUSE_DOUBLE_PRECISION=ON . (plus other variables)
#loop over different settings
for l in range(0, loewe_defs.nt_ns):
	for i in range(0, loewe_defs.cpu_gpu):
		for j in range(0, loewe_defs.double_single):
			for k in range (0, loewe_defs.reconstruct):
				#go into specific folder
				newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
				os.chdir(cur_dir + '/' + newdir)
				#make new files
				os.remove('globaldefs.h')
				f = open('./globaldefs.h', 'w')
				f.write(globaldefs1 + '\n#define NTIME ' + loewe_defs.Nt[l] + '\n#define NSPACE ' + loewe_defs.Ns[l] + '\n' + globaldefs2)
				f.close()
				#make new dir
				newdir2 = 'benchmark_' + loewe_defs.idx[l]
				#cp build dir to new name
				os.system('rm -r -f ' + newdir2)
				#this works only if the new dir does not exist
				shutil.copytree('build', newdir2)
				os.chdir(cur_dir + '/' + newdir + '/' + newdir2 + '/')
				os.system('make clean')
 				os.system('make')
				#send job script
				f = open('./jobscript.sh', 'w')
				f.write(jobscript + loewe_defs.jobscript_switch[i] + jobscript2 + loewe_defs.idx[l])
				f.close()
				os.system('qsub jobscript.sh' );
				os.chdir(cur_dir)
#os.system('mv time_* ./times');

  

