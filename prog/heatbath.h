/** @file
 *
 * Everything required by heatbath's main()
 */
#ifndef _HEATBATHH_
#define _HEATBATHH_
//should only be included in main prog

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_geometry.h"
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_random.h"
#include "host_update_heatbath.h"
#include "host_use_timer.h"
#include "gaugefield_heatbath.h"
#include "opencl_heatbath.h"
#include "logger.hpp"

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

string const version = "0.1b";
using namespace std;

//global random number generator
Random rnd (seed);

//these can be used to measure times on the host
usetimer polytime;
usetimer plaqtime;

//couple of timers
usetimer total_timer;
usetimer init_timer;
usetimer perform_timer;

usetimer copy_to_from_dev_timer;
usetimer copy_on_dev_timer;
usetimer solver_timer;

//to save gaugeobservables
hmc_float plaq, splaq, tplaq;
hmc_complex pol;

void heatbath_time_output(usetimer * total, usetimer * init_timer, usetimer * perform_timer, usetimer * copy_to_from_dev_timer, usetimer* copy_on_dev_timer){

	uint64_t totaltime = (*total).getTime();

	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time
	
	uint64_t init_time = (*init_timer).getTime();
	uint64_t perform_time = (*perform_timer).getTime();
	uint64_t copy2_time = (*copy_on_dev_timer).getTime();
	uint64_t copy1_time = (*copy_to_from_dev_timer).getTime();

	int copy1_steps;
	int copy2_steps;
	
	uint64_t copy1_avgtime;
	uint64_t copy2_avgtime;
	
	copy2_steps = (*copy_on_dev_timer).getNumMeas();
	copy1_steps = (*copy_to_from_dev_timer).getNumMeas();
	
	copy2_avgtime = divide(copy2_time, copy2_steps);
	copy1_avgtime = divide(copy1_time, copy1_steps);

	logger.trace() << "## *******************************************************************";
	logger.trace() << "## Times [mus]:";
	logger.trace() << "## *******************************************************************";
	logger.trace() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
	logger.trace() << "## Total:\t" << setfill(' ') << setw(12) << totaltime;
	logger.trace() << "## Init.:\t" << setfill(' ') << setw(12) << init_time << '\t'<< fixed << setw(5) << setprecision(1) << percent(init_time, totaltime) ;
	logger.trace() << "## Perf.:\t" << setfill(' ') << setw(12) << perform_time << '\t'<< fixed << setw(5) << setprecision(1) << percent(perform_time, totaltime) ;
	logger.trace() << "## *******************************************************************";
	logger.trace() << "## Other:\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg"<< '\t' << setw(5) << "perc";
	logger.trace() << "## Copy1:\t" << setfill(' ') << setw(12) << copy1_time << '\t' << setw(12) << copy1_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(copy1_time, totaltime);
	logger.trace() << "## Copy2:\t" << setfill(' ') << setw(12) << copy2_time << '\t' << setw(12) << copy2_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(copy2_time, totaltime);
	logger.trace() << "## *******************************************************************";

	logger.trace() << "## No output to file implemented yet...";
	/** @todo output to file is not implemented */

	return;
}


#endif /* _HEATBATH_ */
