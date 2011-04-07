#ifndef _USETIMERH_
#define _USETIMERH_

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
//#include <boost/lexical_cast.hpp>

#include "globaldefs.h"
#include "host_input.h"
#include "hmcerrs.h"
#include "klepsydra/klepsydra.hpp"

extern char * benchmark_id;

class usetimer {
  public:
  usetimer() {time_measurement = 0;num_meas = 0;};
  void reset();
  void getTimeAndReset();
  void add();
  void zero();
  uint64_t getTime();
  int getNumMeas();
 private:
  uint64_t time_measurement;
  klepsydra::Monotonic timer;
  int num_meas;
};

void time_output(
	usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy
#ifdef _FERMIONS_
, usetimer * inittimer, usetimer* singletimer, usetimer *Mtimer, usetimer *copytimer, usetimer *scalarprodtimer, usetimer *latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer
#endif
#ifdef _PERFORM_BENCHMARKS_
, int steps
#endif
#ifdef _USEHMC_
, usetimer * hmctimer, usetimer * leapfrogtimer, usetimer * hmcinittimer, usetimer * metropolistimer
#endif
	);


#endif
