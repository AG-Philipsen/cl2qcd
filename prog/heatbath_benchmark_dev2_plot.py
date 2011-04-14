# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil
import loewe_defs


#(benchmark_id) NTIME   NSPACE   VOL4D   totaltime   inittimer   polytime   plaqtime   updatetime   overrelaxtime  (all times average per time-measurement)

s1 = """set term postscript eps color
set key left top
#set size 0.75,0.75
#set xtics 0.05
#set xrange[0.09:0.41]
set xlabel \"V\"
set autoscale y

"""

s2 = """
set output \"Heatbath_Benchmark_Plot_"""
s3 = """.eps\"
set ylabel \""""
s4 = """
set xrange[0:"""
options = """ with linespoints  lw 5 """

cur_dir = os.getcwd()
data_dir = cur_dir + '/heatbath_benchmark/'

cter = 0
#there are at most 8 files
filearray = []
titlearray = []
yaxis = 'Totaltime [ms]', 'Inittime [ms]', 'Copytime [ms]', 'Polyakovloop [ms]', 'Plaquette [ms]', 'Heatbath [ms]', 'Overrelaxation [ms]'

os.chdir(data_dir)

#make plots with all data in it
for i in range(loewe_defs.cpu_gpu_start, loewe_defs.cpu_gpu):
	for j in range(loewe_defs.double_single_start, loewe_defs.double_single):
		for k in range (loewe_defs.reconstruct_start, loewe_defs.reconstruct):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
f = open('./plot', 'w')
f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
for l in range(loewe_defs.obs_start, loewe_defs.obs):
	f.write(s2 + str(l) + s3 + yaxis[l] + '\"\n plot ')
	for m in range(0, cter):
		#divide by 1000 to get to millisecs
		f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
		if(m<cter-1):
			f.write(', ')
f.close()
os.system('gnuplot plot')
os.system('rm plot')

#make plots only with CPU
cter = 0
filearray = []
titlearray = []
if (loewe_defs.cpu_gpu_start == 0):
	i = 0
	for j in range(loewe_defs.double_single_start, loewe_defs.double_single):
		for k in range (loewe_defs.reconstruct_start, loewe_defs.reconstruct):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
	f = open('./plot', 'w')
	f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
	for l in range(loewe_defs.obs_start, loewe_defs.obs):
		f.write(s2 + loewe_defs.cpu_gpu_folder[i] + '_' + str(l) + s3 + yaxis[l] + '\"\n plot ')
		for m in range(0, cter):
			#divide by 1000 to get to millisecs
			f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
			if(m<cter-1):
				f.write(', ')
	f.close()
	os.system('gnuplot plot')
	os.system('rm plot')

#make plots only with GPU
cter = 0
filearray = []
titlearray = []
if (loewe_defs.cpu_gpu == 2):
	i = 1
	for j in range(loewe_defs.double_single_start, loewe_defs.double_single):
		for k in range (loewe_defs.reconstruct_start, loewe_defs.reconstruct):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
	f = open('./plot', 'w')
	f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
	for l in range(loewe_defs.obs_start, loewe_defs.obs):
		f.write(s2 + loewe_defs.cpu_gpu_folder[i] + '_' + str(l) + s3 + yaxis[l] + '\"\n plot ')
		for m in range(0, cter):
			#divide by 1000 to get to millisecs
			f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
			if(m<cter-1):
				f.write(', ')
	f.close()
	os.system('gnuplot plot')
	os.system('rm plot')

#make plots only with double
cter = 0
filearray = []
titlearray = []
if (loewe_defs.double_single_start == 0):
	j = 0
	for i in range(loewe_defs.cpu_gpu_start, loewe_defs.cpu_gpu):
		for k in range (loewe_defs.reconstruct_start, loewe_defs.reconstruct):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
	f = open('./plot', 'w')
	f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
	for l in range(loewe_defs.obs_start, loewe_defs.obs):
		f.write(s2 + loewe_defs.double_single_folder[j] + '_' + str(l) + s3 + yaxis[l] + '\"\n plot ')
		for m in range(0, cter):
			#divide by 1000 to get to millisecs
			f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
			if(m<cter-1):
				f.write(', ')
	f.close()
	os.system('gnuplot plot')
	os.system('rm plot')

#make plots only with single
cter = 0
filearray = []
titlearray = []
if (loewe_defs.double_single == 2):
	j = 1
	for i in range(loewe_defs.cpu_gpu_start, loewe_defs.cpu_gpu):
		for k in range (loewe_defs.reconstruct_start, loewe_defs.reconstruct):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
	f = open('./plot', 'w')
	f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
	for l in range(loewe_defs.obs_start, loewe_defs.obs):
		f.write(s2 + loewe_defs.double_single_folder[j] + '_' + str(l) + s3 + yaxis[l] + '\"\n plot ')
		for m in range(0, cter):
			#divide by 1000 to get to millisecs
			f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
			if(m<cter-1):
				f.write(', ')
	f.close()
	os.system('gnuplot plot')
	os.system('rm plot')

#make plots only with no reconstruct
cter = 0
filearray = []
titlearray = []
if (loewe_defs.reconstruct_start == 0):
	k = 0
	for i in range(loewe_defs.cpu_gpu_start, loewe_defs.cpu_gpu):
		for j in range(loewe_defs.double_single_start, loewe_defs.double_single):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
	f = open('./plot', 'w')
	f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
	for l in range(loewe_defs.obs_start, loewe_defs.obs):
		f.write(s2 + loewe_defs.reconstruct_folder[k] + '_' + str(l) + s3 + yaxis[l] + '\"\n plot ')
		for m in range(0, cter):
			#divide by 1000 to get to millisecs
			f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
			if(m<cter-1):
				f.write(', ')
	f.close()
	os.system('gnuplot plot')
	os.system('rm plot')

#make plots only with reconstruct
cter = 0
filearray = []
titlearray = []
if (loewe_defs.reconstruct == 2):
	k = 1
	for i in range(loewe_defs.cpu_gpu_start, loewe_defs.cpu_gpu):
		for j in range(loewe_defs.double_single_start, loewe_defs.double_single):
			cter+=1
			newdir = loewe_defs.cpu_gpu_folder[i] + '_' + loewe_defs.double_single_folder[j] + '_' + loewe_defs.reconstruct_folder[k]
			filename = 'time_B_' + newdir
			title = '\"' + loewe_defs.cpu_gpu_folder[i] + loewe_defs.double_single_folder[j] + loewe_defs.reconstruct_folder[k]+ '\"'
			filearray.append(filename)
			titlearray.append(title)
	f = open('./plot', 'w')
	f.write(s1 + s4 + loewe_defs.Vol[loewe_defs.nt_ns - 1] + ']\n')
	for l in range(loewe_defs.obs_start, loewe_defs.obs):
		f.write(s2 + loewe_defs.reconstruct_folder[k] + '_' + str(l) + s3 + yaxis[l] + '\"\n plot ')
		for m in range(0, cter):
			#divide by 1000 to get to millisecs
			f.write('\"./' + filearray[m] + '\" using 5:($' + str(l+1+5) + '/1000.) ' + options + ' title ' + titlearray[m])
			if(m<cter-1):
				f.write(', ')
	f.close()
	os.system('gnuplot plot')
	os.system('rm plot')

