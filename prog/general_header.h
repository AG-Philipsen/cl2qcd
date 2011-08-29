/** @file
 *
 * Header files that are included in every executable.
 */
#ifndef _GENERALHEADERH_
#define _GENERALHEADERH_

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_use_timer.h"
#include "host_random.h"

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#include "logger.hpp"

using namespace std;

//global random number generator
Random rnd (seed);

//these can be used to measure times on the host spent on gaugeobservables
usetimer poly_timer;
usetimer plaq_timer;

//timers for overall times
usetimer total_timer;
usetimer init_timer;
usetimer perform_timer;

//to save gaugeobservables on host
hmc_float plaq, splaq, tplaq;
hmc_complex pol;

//output-method for timer above and the copy-timer from (each) device
///@todo For hybrid case one will have to add more timer here
void general_time_output(usetimer * total, usetimer * init_timer, usetimer * perform_timer, usetimer * copy_to_from_dev_timer, usetimer* copy_on_dev_timer, usetimer * plaq_timer, usetimer * poly_timer){

	uint64_t totaltime = (*total).getTime();

	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time
	
	uint64_t init_time = (*init_timer).getTime();
	uint64_t perform_time = (*perform_timer).getTime();
	uint64_t copy2_time = (*copy_on_dev_timer).getTime();
	uint64_t copy1_time = (*copy_to_from_dev_timer).getTime();
	uint64_t plaq_time = (*plaq_timer).getTime();
	uint64_t poly_time = (*poly_timer).getTime();

	int copy1_steps = (*copy_to_from_dev_timer).getNumMeas();
	int copy2_steps = (*copy_on_dev_timer).getNumMeas();
	int plaq_steps = (*plaq_timer).getNumMeas();
	int poly_steps = (*poly_timer).getNumMeas();
	
	uint64_t copy1_avgtime = divide(copy1_time, copy1_steps);
	uint64_t copy2_avgtime = divide(copy2_time, copy2_steps);
	uint64_t poly_avgtime = divide(poly_time, poly_steps);
	uint64_t plaq_avgtime = divide(plaq_time, plaq_steps);

	logger.trace() << "## *******************************************************************";
	logger.trace() << "## General Times [mus]:";
	logger.trace() << "## *******************************************************************";
	logger.trace() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
	logger.trace() << "## Total:\t" << setfill(' ') << setw(12) << totaltime;
	logger.trace() << "## Init.:\t" << setfill(' ') << setw(12) << init_time << '\t'<< fixed << setw(5) << setprecision(1) << percent(init_time, totaltime) ;
	logger.trace() << "## Perf.:\t" << setfill(' ') << setw(12) << perform_time << '\t'<< fixed << setw(5) << setprecision(1) << percent(perform_time, totaltime) ;
	logger.trace() << "## *******************************************************************";
	logger.trace() << "## Other:\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg"<< '\t' << setw(5) << "perc";
	logger.trace() << "## Copy1:\t" << setfill(' ') << setw(12) << copy1_time << '\t' << setw(12) << copy1_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(copy1_time, totaltime);
	logger.trace() << "## Copy2:\t" << setfill(' ') << setw(12) << copy2_time << '\t' << setw(12) << copy2_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(copy2_time, totaltime);
	logger.trace() << "## Plaq.:\t" << setfill(' ') << setw(12) << plaq_time << '\t' << setw(12) << plaq_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(plaq_time, totaltime);
	logger.trace() << "## Poly.:\t" << setfill(' ') << setw(12) << poly_time << '\t' << setw(12) << poly_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(poly_time, totaltime);
	logger.trace() << "## *******************************************************************";

	logger.trace() << "## No output of times to file implemented yet...";
	/** @todo output to file is not implemented */
	//See older files for example code

	return;
}

#endif /* _<GENERALHEADERH_ */
