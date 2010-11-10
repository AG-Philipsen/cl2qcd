PROGRAM = hmc
PROGNAME = hmc

GLOBALS=globals.h types.h hmcerrs.h testing.h

CXX = g++
CXXFLAGS = -Wall -pedantic -D_RECONSTRUCT_TWELVE_
#CXXFLAGS = -Wall -pedantic

LDFLAGS =
LIBS = -lm 

INCLUDES =

OBJECTS = operations.o gaugeobservables.o testing.o geometry.o input.o

main: $(PROGRAM).cpp $(GLOBALS) $(OBJECTS)
	$(CXX) -o ${PROGNAME} $(CXXFLAGS) $(PROGRAM).cpp ${OBJECTS} $(LDFLAGS) $(LIBS)

%.o: %.cpp %.h $(GLOBALS)
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f hmc *.o

.PHONY: all clean
