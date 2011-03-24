# -*- coding: utf-8 -*-

#!/usr/bin/env python

#variables for loops
cpu_gpu = 2
double_single = 2
reconstruct = 2
nt_ns = 4

#arrays with values for the Makefile and globaldefs.h
Nt = '4', '4', '4', '4', '4',  '4',  '4',  '4',  '6',  '6',  '6',  '8',  '8',  '8',  '8',  '8',  '8',  '8',  '10',  '10',  '12',  '14'
Ns = '8', '12', '14', '16', '18',  '20',  '22',  '24',  '12',  '24',  '36',  '16',  '18',  '20',  '22',  '24',  '32',  '42',  '20',  '40',  '24',  '28'
cpu_gpu_text = '#DEFGPU=-D_USEGPU_\n', 'DEFGPU=-D_USEGPU_\n'
double_single_text = '#DEFDOUBLE=-D_USEDOUBLEPREC_\n', 'DEFDOUBLE=-D_USEDOUBLEPREC_\n'
reconstruct_text = '#DEFREC=-D_RECONSTRUCT_TWELVE_\n', 'DEFREC=-D_RECONSTRUCT_TWELVE_\n'
idx = '001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019', '020', '021', '022', 
cpu_gpu_folder = 'C', 'G'
double_single_folder = 'D', 'S'
reconstruct_folder = 'N', 'R'
jobscript_switch = 'default', 'gpu'