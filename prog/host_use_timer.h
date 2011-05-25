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
#include "inputparameters.h"
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
	 * Retrieve the aggregated measured time in microseconds (10^6s).
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

/* 
LZ: old version, relies on compiler flags which will be removed; keep copy for now so that we will eventually be able to copy and paste parts of the code in the function time_output_hmc which has to be written later

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
*/

/**
 * Print statistics given by the passed timers to stdout and to a file
 * Designed for use in heatbath.cpp
 */
void time_output_heatbath(usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy);
/**
 * Print statistics given by the passed timers to stdout and to a file
 * Designed for use in inverter.cpp
 */
void time_output_inverter(usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy, usetimer * inittimer, usetimer* singletimer, usetimer *Mtimer, usetimer *copytimer, usetimer *scalarprodtimer, usetimer *latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer);

#endif
