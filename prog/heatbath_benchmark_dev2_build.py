# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil
import dev2_defs


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
#define GAUGEMOMENTASIZE NDIM*VOL4D*(NC*NC-1)
#define GAUGEFIELDSIZE NC*NC*NDIM*VOL4D

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

#suppose you call the script from the folder /path-to-hmc/
#	and the prog is lying in /path-to-hmc/prog
cur_dir = os.getcwd()
prog_dir = cur_dir + '/prog/'
#make dir for data-collection later
data_dir = cur_dir + '/heatbath_benchmark/'
#check if path exists
if not os.path.exists(data_dir):
	os.makedirs(data_dir)

#collect options for cmake
cmakeopt = '  -DPERFORM_BENCHMARKS=ON ' #perform benchmarks
cmakeopt += ' -DCMAKE_BUILD_TYPE=Release ' #"Release" setting
cmakeopt += ' -DCMAKE_MODULE_PATH=/home/pinke/.cmake/modules' #path to find FindOpenCL.cmake

#this script builds a hmc for each setting
#go to prog dir
os.chdir(prog_dir)
#loop over different settings
for i in range(dev2_defs.cpu_gpu_start, dev2_defs.cpu_gpu):
	for j in range(dev2_defs.double_single_start, dev2_defs.double_single):
		for k in range (dev2_defs.reconstruct_start, dev2_defs.reconstruct):
			#collect cmake options
			cmakeopt += dev2_defs.cpu_gpu_opt[i]
			cmakeopt += dev2_defs.double_single_opt[j]
			cmakeopt += dev2_defs.reconstruct_opt[k]
			for l in range(dev2_defs.nt_ns_start, dev2_defs.nt_ns):
				#make new globaldef.h file
				#os.remove('globaldefs.h')
				f = open('./globaldefs.h', 'w')
				f.write(globaldefs1 + '\n#define NTIME ' + dev2_defs.Nt[l] + '\n#define NSPACE ' + dev2_defs.Ns[l] + '\n' + globaldefs2)
				f.close()
				#prepare folder for program
				newdir = dev2_defs.cpu_gpu_folder[i] + '_' + dev2_defs.double_single_folder[j] + '_' + dev2_defs.reconstruct_folder[k] + '_' + dev2_defs.idx[l]
				#check if path exists
				if not os.path.exists(newdir):
					os.makedirs(newdir)
				#go into specific folder
				os.chdir(prog_dir + '/' + newdir)
				#prepare makefiel via cmake
				os.system('cmake ' + cmakeopt + ' ..')
				os.system('make')
				os.chdir(prog_dir)
