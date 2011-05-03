/** @file
 * Time measurement utilities
 */
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

/**
 * A wrapper for the monotonic Klepsydra timer
 * that can add up the time of multiple measurements.
 *
 * @todo Make this a child of klepsydra::Monotonic
 * @todo Merge this functionality back into Klepsydra.
 */
class usetimer {
  public:
  /**
   * Default and only constructor.
   */
  usetimer() : time_measurement(0), num_meas(0) { };
  /**
   * Reset the current time measurement.
   */
  void reset();
  /**
   * Store passed time and reset.
   */
  void getTimeAndReset();
  /**
   * Add passed time to any previously measured time and reset.
   */
  void add();
  /**
   * Reset the aggregated measurement information.
   */
  void zero();
  /**
   * Retrieve the aggregated measured time.
   */
  uint64_t getTime();
  /**
   * Retrieve the number of measurements performed.
   */
  int getNumMeas();
 private:
  /**
   * The aggregated measured time.
   */
  uint64_t time_measurement;
  /**
   * Currrently running measurement (not included in time_measurement)
   */
  klepsydra::Monotonic timer;
  /**
   * The number of measurements aggregated so far
   */
  int num_meas;
};

/**
 * Print statistics given by the passed timers to stdout
 */
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


void time_heatbath(usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy);


#endif
