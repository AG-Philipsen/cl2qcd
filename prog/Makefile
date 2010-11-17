PROGRAM = hmc
PROGNAME = hmc

GLOBALS=globals.h types.h hmcerrs.h testing.h

#CP:
#to supress deprecated warnings in the tmlqcd read functions: -Wno-write-strings
CXX = g++
CXXFLAGS = -Wall -pedantic $(xml2-config --cflags) $(xml2-config --libs) -Wno-write-strings
#CXXFLAGS = -Wall -pedantic -Wno-write-strings -D_RECONSTRUCT_TWELVE_

#LDFLAGS = -L/home/pinke/lime/lime/lib
LDFLAGS = -L/home/zeidlewicz/bin/lime/lime-lib/lib
LIBS = -lm -llime -lxml2

#INCLUDES = -I/home/pinke/lime/lime/include -I/usr/include/libxml2
INCLUDES = -I/home/zeidlewicz/bin/lime/lime-lib/include -I/usr/include/libxml2

OBJECTS = operations.o gaugeobservables.o testing.o geometry.o input.o readgauge.o

main: $(PROGRAM).cpp $(GLOBALS) $(OBJECTS) 
	$(CXX) -o ${PROGNAME} $(CXXFLAGS) $(PROGRAM).cpp ${OBJECTS} $(LDFLAGS) $(LIBS) $(INCLUDES) 

#readgauge.o: readgauge.cpp readgauge.h readtmlqcd.o $(GLOBALS)
#	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) readgauge.cpp

%.o: %.cpp %.h $(GLOBALS)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBS) $<

clean:
	rm -f hmc *.o

.PHONY: all clean
