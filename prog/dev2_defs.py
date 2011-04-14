# -*- coding: utf-8 -*-

#!/usr/bin/env python

#variables for loops
cpu_gpu = 2
cpu_gpu_start = 1
double_single = 2
double_single_start = 1
reconstruct = 2
reconstruct_start = 1
nt_ns = 1
nt_ns_start = 0
obs_start = 0
obs = 7

#arrays with values for the Makefile and globaldefs.h
Nt = '4', '4', '4', '4', '4',  '4',  '4',  '4',  '6',  '6',  '6',  '8',  '8',  '8',  '8',  '8',  '8',  '8',  '10',  '10',  '12',  '14'
Ns = '8', '12', '14', '16', '18',  '20',  '22',  '24',  '12',  '24',  '36',  '16',  '18',  '20',  '22',  '24',  '32',  '42',  '20',  '40',  '24',  '28'
#Vol = '2048', '6912', '10976', '16384', '23328', '32000', '42592', '55296', '10368', '82944', '279936', '32768', '46656', '64000', '85184', '110592', '262144', '592704', '80000', '640000', '165888', '307328'
Vol = '2200', '7200', '12000', '17000', '24000', '34000', '44000', '56000', '56000', '84000', '280000', '280000','280000','280000','280000','280000','280000', '600000', '600000', '650000', '650000', '650000'
cpu_gpu_text = '#DEFGPU=-D_USEGPU_\n', 'DEFGPU=-D_USEGPU_\n'
double_single_text = '#DEFDOUBLE=-D_USEDOUBLEPREC_\n', 'DEFDOUBLE=-D_USEDOUBLEPREC_\n'
reconstruct_text = '#DEFREC=-D_RECONSTRUCT_TWELVE_\n', 'DEFREC=-D_RECONSTRUCT_TWELVE_\n'
idx = '001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019', '020', '021', '022', 
cpu_gpu_folder = 'C', 'G'
double_single_folder = 'D', 'S'
reconstruct_folder = 'N', 'R'
jobscript_switch = 'default', 'gpu'
cpu_gpu_opt = ' -DUSE_GPU=OFF', ' -DUSE_GPU=ON'
double_single_opt = ' -DUSE_DOUBLE_PRECISION=ON', ' -DUSE_DOUBLE_PRECISION=OFF'
reconstruct_opt = ' -DUSE_RECONSTRUCT_TWELVE=OFF', ' -DUSE_RECONSTRUCT_TWELVE=ON'