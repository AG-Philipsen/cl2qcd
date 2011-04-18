# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil
import dev2_defs

#suppose you call the script from the folder /path-to-hmc/
#	and the prog is lying in /path-to-hmc/prog
cur_dir = os.getcwd()
prog_dir = cur_dir + '/prog/'
data_dir = cur_dir + '/heatbath_benchmark/'

os.chdir(data_dir)
data_dir += '/' + str(dev2_defs.benchmark_number)
#make sure data folder exists
if not os.path.exists(data_dir):
	os.makedirs(data_dir)

#loop over different settings
for i in range(dev2_defs.cpu_gpu_start, dev2_defs.cpu_gpu):
	for j in range(dev2_defs.double_single_start, dev2_defs.double_single):
		for k in range (dev2_defs.reconstruct_start, dev2_defs.reconstruct):
			newdir = dev2_defs.cpu_gpu_folder[i] + '_' + dev2_defs.double_single_folder[j] + '_' + dev2_defs.reconstruct_folder[k]
			#go into data folder
			os.chdir(data_dir)
			#make empty data-files
			filename = 'time_B_' + newdir
			os.system('rm -f ' + filename)
			f = open('./' + filename, 'w')
			f.close()
			#copy time_output into data_dir
			os.chdir(prog_dir)
			for l in range(dev2_defs.nt_ns_start, dev2_defs.nt_ns):
				#go into specific folder
				newdir_spec = newdir + '_' + dev2_defs.idx[l]
				filename_spec = 'time_B_' + newdir_spec
				os.chdir(newdir_spec)
				#cp file to data_dir
				os.system('cp ' + filename_spec + ' ' + data_dir) 
				os.chdir(prog_dir)
			#collect data
			os.chdir(data_dir)
			for l in range(dev2_defs.nt_ns_start, dev2_defs.nt_ns):
				newdir_spec = newdir + '_' + dev2_defs.idx[l]
				filename_spec = 'time_B_' + newdir_spec
				newdir2 = 'benchmark_' + dev2_defs.idx[l]
				#os.chdir(cur_dir + '/' + newdir + '/' + newdir2 + '/')
				os.system('cat ' + filename_spec + ' >>  ' + filename) 
				#remove specific data file
				os.system('rm ' + filename_spec)

  

