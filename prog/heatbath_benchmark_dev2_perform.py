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
#make dir for data-collection later
data_dir = cur_dir + '/heatbath_benchmark/'

#this runs the hmc for different setting. the executable is supposed to be built beforehand
#go to prog dir
os.chdir(prog_dir)
#loop over different settings
for i in range(dev2_defs.cpu_gpu_start, dev2_defs.cpu_gpu):
	for j in range(dev2_defs.double_single_start, dev2_defs.double_single):
		for k in range (dev2_defs.reconstruct_start, dev2_defs.reconstruct):
			for l in range(dev2_defs.nt_ns_start, dev2_defs.nt_ns):
				#go in folder for program
				newdir = dev2_defs.cpu_gpu_folder[i] + '_' + dev2_defs.double_single_folder[j] + '_' + dev2_defs.reconstruct_folder[k] + '_' + dev2_defs.idx[l]
				##go into specific folder
				os.chdir(prog_dir + '/' + newdir)
				os.system('./hmc ' + cur_dir + '/input_heatbath_benchmark ' + str(dev2_defs.idx[l]))
				os.chdir(prog_dir)