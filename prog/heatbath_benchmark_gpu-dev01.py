# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil


#test hmc with _D_RECONSTRUCT_TWELVE

#here, only settings for gpu-dev01 are included!!
makefile = """

PROGRAM = hmc
VERSION=0.1
PROGNAME = hmc.v$(VERSION)

CXX = g++
#gpu-dev01:
LIMEPATH=/home/pinke/lime/lime/
XMLPATH=/home/pinke/libxml2/libxml2
ATISDKPATH=/home/pinke/ati-stream-sdk-v2.3-lnx64

LIMEINCDIR=$(LIMEPATH)/include 
LIMELIBDIR=$(LIMEPATH)/lib
ATISDKINCDIR=$(ATISDKPATH)/include
ATISDKLIBDIR=$(ATISDKPATH)/lib/x86_64/
XMLINCDIR=$(XMLPATH)
#gpudev01:
XMLINCDIR=$(XMLPATH)/include/libxml2
XMLLIBDIR=$(XMLPATH)/lib
DEFS=$(DEFGPU) $(DEFDOUBLE) $(DEFREC) $(DEFOMP) -DSOURCEDIR=$(OCLSOURCEDIR)
INCLUDES = -I$(LIMEINCDIR) -I$(XMLINCDIR) -I$(ATISDKINCDIR) -I$(HMCPATH)/host_operations
LDFLAGS = -L$(LIMELIBDIR) -L$(ATISDKLIBDIR) -L$(XMLLIBDIR)
LIBS = -lm -llime -lxml2 -lOpenCL -fopenmp

CXXFLAGS = -Wall -pedantic -std=gnu++0x $(DEFS)

GLOBALS=globaldefs.h types.h hmcerrs.h

OBJECTS = host_operations_complex.o \
	host_operations_gaugefield.o \
	host_operations_spinor.o \
	host_gaugeobservables.o \
	host_testing.o \
	host_geometry.o \
	host_input.o \
	host_readgauge.o \
	host_random.o \
	host_update_heatbath.o \
	host_timer.o \
	opencl.o \
	host_use_timer.o \
	host_gaugefieldoperations.o \
	host_writegaugefield.o \
	host_solver.o \
	host_fermionobservables.o


main: $(PROGRAM).cpp $(GLOBALS) $(PROGRAM).h $(OBJECTS)
	$(CXX) -o ${PROGNAME} $(CXXFLAGS) $(PROGRAM).cpp ${OBJECTS} $(LDFLAGS) $(LIBS) $(INCLUDES) 

host_gaugefieldoperations.o: host_gaugefieldoperations.cpp host_gaugefieldoperations.h $(GLOBALS) host_use_timer.o  host_readgauge.o host_input.o host_gaugeobservables.o host_operations_complex.o
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

opencl.o: opencl.cpp opencl.h $(GLOBALS)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

use_timer.o: host_use_timer.cpp host_use_timer.h $(GLOBALS) host_timer.o
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $<

%.o: %.cpp %.h $(GLOBALS)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES)  $<

clean:
	rm -f $(PROGNAME) *.o

.PHONY: all clean
"""

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
#define NUMTHREADS 1
#else
#define NUMTHREADS 1
#endif

#endif
"""

#variables for loops
cpu_gpu = 1
double_single = 1
reconstruct = 1
nt_ns = 1

#arrays with values for the Makefile and globaldefs.h
Nt = '4', '4', '4', '4', '4',  '4',  '4',  '4',  '6',  '6',  '6',  '8',  '8',  '8',  '8',  '8',  '8',  '8',  '10',  '10',  '12',  '14'
Ns = '8', '12', '14', '16', '18',  '20',  '22',  '24',  '12',  '24',  '36',  '16',  '18',  '20',  '22',  '24',  '32',  '42',  '20',  '40',  '24',  '28'
cpu_gpu_text = '#DEFGPU=-D_USEGPU_\n', 'DEFGPU=-D_USEGPU_\n'
double_single_text = '#DEFDOUBLE=-D_USEDOUBLEPREC_\n', 'DEFDOUBLE=-D_USEDOUBLEPREC_\n'
reconstruct_text = '#DEFREC=-D_RECONSTRUCT_TWELVE_\n', 'DEFREC=-D_RECONSTRUCT_TWELVE_\n'
idx = '001 ', '002 ', '003 ', '004 ', '005 ', '006 ', '007 ', '008 ', '009 ', '010 ', '011 ', '012 ', '013 ', '014 ', '015 ', '016 ', '017 ', '018 ', '019 ', '020 ', '021 ', '022 ', 

#compile

cur_dir = os.getcwd()

#loop over different settings
for l in range(0, nt_ns):
	for i in range(0, cpu_gpu):
		for j in range(0, double_single):
			for k in range (0, reconstruct):
				#make the two files
				os.remove('Makefile')
				f = open('./Makefile', 'w')
				f.write('OCLSOURCEDIR=\"\\\"/home/pinke/hmc/benchmark_' + str(i) + '_' + str(j) + '_' + str(l) + '_' + idx[l] + '\"\\\"\n')
				f.write(cpu_gpu_text[i] + double_single_text[j] + reconstruct_text[k] + makefile)
				f.close()
				os.remove('globaldefs.h')
				f = open('./globaldefs.h', 'w')
				f.write(globaldefs1 + '\n#define NTIME ' + Nt[l] + '\n#define NSPACE ' + Ns[l] + '\n' + globaldefs2)
				f.close()
				os.system('make clean')
 				os.system('make')
				#make new dir
				newdir = 'benchmark_' + str(i) + '_' + str(j) + '_' + str(l) + '_' + idx[l]
				os.system('rm -r ' + newdir)
				os.system('mkdir ' + newdir)
				os.system('cp *.cl ./' + newdir)
				os.system('cp types.h ./' + newdir)
				os.system('cp globaldefs.h ./' + newdir)
				os.chdir(cur_dir + '/' + newdir)
				#execute hmc
				os.system('./hmc.v0.1 input_heatbath_benchmark ' + idx[l] );
				os.chdir(cur_dir)
#os.system('mv time_* ./times');

  

