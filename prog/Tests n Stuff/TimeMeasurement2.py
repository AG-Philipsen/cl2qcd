# -*- coding: utf-8 -*-

#!/usr/bin/env python
import os
import sys
import shutil


#test hmc with _D_RECONSTRUCT_TWELVE

makefile = """
PROGRAM = hmc
PROGNAME = hmc

CXX = g++

DEFS=-D_RECONSTRUCT_TWELVE_
#DEFS=-D_RECONSTRUCT_TWELVE_ -D_OPENMP

INCLUDES = -I/home/pinke/lime/lime/include -I/usr/include/libxml2
#INCLUDES = -I/home/zeidlewicz/bin/lime/lime-lib/include -I/usr/include/libxml2 -I/home/zeidlewicz/bin/ati-stream-sdk-v2.2-lnx64/include/

LDFLAGS = -L/home/pinke/lime/lime/lib
#LDFLAGS = -L/home/zeidlewicz/bin/lime/lime-lib/lib -I/home/zeidlewicz/bin/ati-stream-sdk-v2.2-lnx64/lib/x86_64/

LIBS = -lm -llime -lxml2 -lOpenCL

CXXFLAGS = -Wall -pedantic $(xml2-config --cflags) $(xml2-config --libs) -std=gnu++0x -fopenmp $(DEFS)

GLOBALS=globaldefs.h types.h hmcerrs.h

OBJECTS = operations.o \
	gaugeobservables.o \
	testing.o \
	geometry.o \
	input.o \
	readgauge.o \
	random.o \
	update.o \
	timer.o \
	opencl.o

main: $(PROGRAM).cpp $(GLOBALS) $(PROGRAM).h $(OBJECTS)
	$(CXX) -o ${PROGNAME} $(CXXFLAGS) $(PROGRAM).cpp ${OBJECTS} $(LDFLAGS) $(LIBS) $(INCLUDES) 

%.o: %.cpp %.h $(GLOBALS)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES)  $<

clean:
	rm -f $(PROGNAME) *.o

.PHONY: all clean
"""

makefile2 = """
PROGRAM = hmc
PROGNAME = hmc

CXX = g++

DEFS=-D_RECONSTRUCT_TWELVE_
#DEFS=-D_RECONSTRUCT_TWELVE_ -D_OPENMP

INCLUDES = -I/home/pinke/lime/lime/include -I/usr/include/libxml2
#INCLUDES = -I/home/zeidlewicz/bin/lime/lime-lib/include -I/usr/include/libxml2 -I/home/zeidlewicz/bin/ati-stream-sdk-v2.2-lnx64/include/

LDFLAGS = -L/home/pinke/lime/lime/lib
#LDFLAGS = -L/home/zeidlewicz/bin/lime/lime-lib/lib -I/home/zeidlewicz/bin/ati-stream-sdk-v2.2-lnx64/lib/x86_64/

LIBS = -lm -llime -lxml2 -lrt

CXXFLAGS = -Wall -pedantic $(xml2-config --cflags) $(xml2-config --libs) -std=gnu++0x -fopenmp $(DEFS)

GLOBALS=globaldefs.h types.h hmcerrs.h

OBJECTS = operations.o \
	gaugeobservables.o \
	testing.o \
	geometry.o \
	input.o \
	readgauge.o \
	random.o \
	update.o \
	timer.o \
	opencl.o

main: $(PROGRAM).cpp $(GLOBALS) $(PROGRAM).h $(OBJECTS)
	$(CXX) -o ${PROGNAME} $(CXXFLAGS) $(PROGRAM).cpp ${OBJECTS} $(LDFLAGS) $(LIBS) $(INCLUDES) 

%.o: %.cpp %.h $(GLOBALS)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES)  $<

clean:
	rm -f $(PROGNAME) *.o

.PHONY: all clean
"""

res = 2
res2 = 6

#compile

txt1 = """#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

//it could be a good idea to define Nt and Ns at compile time
//usually you stick to one volume for quite a while anyways...
#define NSPACE """

txt2 = """
#define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

#define su2_entries 4

#define PI 	3.14159265358979

#endif
"""
bla = os.getcwd()

#make new makefile

os.remove('Makefile')
f = open('./Makefile', 'w')
f.write(makefile)
f.close()

#executables must exist!
for i in range (0, res2):
 os.system('make clean')
 os.remove('./globaldefs.h')
 f = open('./globaldefs.h', 'w')
 f.write(txt1)
 f.write(str((i+1)*4) + '\n#define NTIME ' + str((i+1)*4))
 f.write(txt2)
 f.close()
 os.system('make')
 os.rename('hmc', 'hmc_' + str((i+1)*4))
 os.chdir(bla + '/TimeTest2')
 os.remove('hmc_' + str((i+1)*4))
 os.chdir(bla)
 shutil.move('hmc_' + str((i+1)*4), bla + '/TimeTest2')
 
#make makefile without _D_RECONSTRUCT_TWELVE

os.remove('Makefile')
f = open('./Makefile', 'w')
f.write(makefile2)
f.close()


#execute
#input files must already exist!

os.chdir(bla + '/TimeTest2')
for i in range (0,res):
 for j in range (0, res2):
  os.system('asdfsdjkf')
  os.system('./hmc_' + str((j+1)*4) + ' input' + str(i+1));
