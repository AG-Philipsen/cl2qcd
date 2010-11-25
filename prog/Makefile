PROGRAM = hmc
PROGNAME = hmc

CXX = g++

DEFS=-D_RECONSTRUCT_TWELVE_
#DEFS=-D_RECONSTRUCT_TWELVE_ -D_OPENMP

#INCLUDES = -I/home/pinke/lime/lime/include -I/usr/include/libxml2
INCLUDES = -I/home/zeidlewicz/bin/lime/lime-lib/include -I/usr/include/libxml2 -I/home/zeidlewicz/bin/ati-stream-sdk-v2.2-lnx64/include/

#LDFLAGS = -L/home/pinke/lime/lime/lib
LDFLAGS = -L/home/zeidlewicz/bin/lime/lime-lib/lib -I/home/zeidlewicz/bin/ati-stream-sdk-v2.2-lnx64/lib/x86_64/

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
