# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil
import loewe_defs

cur_dir = os.getcwd()

#suppose there is a dir for every option-setting (8 dir) with a suiting build dir in it
#loop over different settings
for i in range(loewe_defs.cpu_gpu_start, loewe_defs.cpu_gpu):
	for j in range(loewe_defs.double_single_start, loewe_defs.double_single):
		for k in range (loewe_defs.reconstruct_start, loewe_defs.reconstruct):
			#go into specific folder
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			os.chdir(cur_dir + '/' + newdir)
			# 'B' stands for benchmarking
			filename = 'time_B_' + newdir
			os.system('rm -f ' + filename)
			f = open('./' + filename, 'w')
			f.close()
			for l in range(loewe_defs.nt_ns_start, loewe_defs.nt_ns):
				newdir2 = 'benchmark_' + loewe_defs.idx[l]
				#os.chdir(cur_dir + '/' + newdir + '/' + newdir2 + '/')
				os.system('cat ' + newdir2 + '/' + filename + '_' + str(loewe_defs.idx[l]) + ' >>  ' + filename) 
			os.system('mv ' + filename + ' ../data');

  

