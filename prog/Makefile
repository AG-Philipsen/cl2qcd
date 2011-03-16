PROGRAM = hmc
VERSION=0.1
PROGNAME = hmc.v$(VERSION)

CXX = g++
XMLPATH=/usr/include/libxml2
HMCPATH=/home/pinke/hmc/prog
LIMEPATH=/home/pinke/lime/lime
OCLSOURCEDIR="\"/home/pinke/hmc/prog"\"
ATISDKPATH=/home/pinke/PHD/OpenCL/ati-stream-sdk-v2.3-lnx64
#OCLSOURCEDIR="\"/home/zeidlewicz/hmc/hmc/prog\""
#LIMEPATH=/home/zeidlewicz/bin/lime/lime-lib
#ATISDKPATH=/home/zeidlewicz/bin/ati-stream-sdk-v2.3-lnx64

#gpu-dev01:
#OCLSOURCEDIR="\"/home/pinke/hmc/prog"\"
#LIMEPATH=/home/pinke/lime/lime/
#XMLPATH=/home/pinke/libxml2/libxml2
#ATISDKPATH=/home/pinke/ati-stream-sdk-v2.3-lnx64

#LOEWE-CSC
#OCLSOURCEDIR="\"/home/hfftheo/zeidlewicz/heatbath-benchmark"\"
#LIMEPATH=/home/hfftheo/zeidlewicz/libs/lime
#XMLPATH=/home/hfftheo/zeidlewicz/libs/libxml2
##ATISDKPATH=/home/hfftheo/zeidlewicz/opencl/ati-stream-sdk-v2.2-lnx64
## this path is not really needed, just for comparability...
#ATISDKPATH=/cm/shared/apps/atistreamsdk/2.3


#DEFGPU=-D_USEGPU_
DEFDOUBLE=-D_USEDOUBLEPREC_
#DEFREC=-D_RECONSTRUCT_TWELVE_
#DEFOMP=-D_OPENMP
#DEFTESTING=-D_TESTING_
DEFFERMIONS=-D_FERMIONS_

LIMEINCDIR=$(LIMEPATH)/include 
LIMELIBDIR=$(LIMEPATH)/lib
ATISDKINCDIR=$(ATISDKPATH)/include
ATISDKLIBDIR=$(ATISDKPATH)/lib/x86_64/
XMLINCDIR=$(XMLPATH)
#gpudev01 OR LOEWE-CSC:
#XMLINCDIR=$(XMLPATH)/include/libxml2
XMLLIBDIR=$(XMLPATH)/lib
DEFS=$(DEFGPU) $(DEFDOUBLE) $(DEFFERMIONS) $(DEFTESTING) $(DEFREC) $(DEFOMP) -DSOURCEDIR=$(OCLSOURCEDIR)
INCLUDES = -I$(LIMEINCDIR) -I$(XMLINCDIR) -I$(ATISDKINCDIR) -I$(HMCPATH)/host_operations
LDFLAGS = -L$(LIMELIBDIR) -L$(ATISDKLIBDIR) -L$(XMLLIBDIR)
LIBS = -lm -llime -lxml2 -lOpenCL -fopenmp

CXXFLAGS = -Wall -pedantic -std=gnu++0x $(DEFS)

GLOBALS=globaldefs.h types.h hmcerrs.h

OBJECTS = host_operations_complex.o \
	host_operations_gaugefield.o \
	host_operations_spinor.o \
	host_operations_spinorfield.o \
	host_operations_fermionmatrix.o \
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
