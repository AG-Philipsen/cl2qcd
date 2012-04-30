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
#include "types.h"
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_use_timer.h"
#include "host_random.h"

#include "exceptions.h"

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
void general_time_output(usetimer * total, usetimer * init_timer, usetimer * perform_timer, usetimer * plaq_timer, usetimer * poly_timer)
{

	uint64_t totaltime = (*total).getTime();

	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time

	uint64_t init_time = (*init_timer).getTime();
	uint64_t perform_time = (*perform_timer).getTime();
	uint64_t plaq_time = (*plaq_timer).getTime();
	uint64_t poly_time = (*poly_timer).getTime();

	int plaq_steps = (*plaq_timer).getNumMeas();
	int poly_steps = (*poly_timer).getNumMeas();

	uint64_t poly_avgtime = divide(poly_time, poly_steps);
	uint64_t plaq_avgtime = divide(plaq_time, plaq_steps);

	logger.info() << "## *******************************************************************";
	logger.info() << "## General Times [mus]:";
	logger.info() << "## *******************************************************************";
	logger.info() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
	logger.info() << "## Total:\t" << setfill(' ') << setw(12) << totaltime;
	logger.info() << "## Init.:\t" << setfill(' ') << setw(12) << init_time << '\t' << fixed << setw(5) << setprecision(1) << percent(init_time, totaltime) ;
	logger.info() << "## Perf.:\t" << setfill(' ') << setw(12) << perform_time << '\t' << fixed << setw(5) << setprecision(1) << percent(perform_time, totaltime) ;
	logger.info() << "## *******************************************************************";
	logger.info() << "## Host-Obs:\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg" << '\t' << setw(5) << "perc";
	logger.info() << "## Plaq.:\t" << setfill(' ') << setw(12) << plaq_time << '\t' << setw(12) << plaq_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(plaq_time, totaltime);
	logger.info() << "## Poly.:\t" << setfill(' ') << setw(12) << poly_time << '\t' << setw(12) << poly_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(poly_time, totaltime);
	logger.info() << "## *******************************************************************";

	logger.info() << "## writing general times to file: \"general_time_output\"";
	ofstream ofile;
	ofile.open("general_time_output");
	if(ofile.is_open()) {
		ofile  << "## *******************************************************************" << endl;
		ofile  << "## General Times [mus]:" << endl;
		ofile << "## Total\tInit\tPerformance\tHost-Plaq\tHost-Pol" << endl;
		ofile  << totaltime << "\t" << init_time << '\t' << perform_time << '\t' << plaq_time << '\t' << poly_time << endl;
		ofile.close();
	} else {
		logger.warn() << "Could not open output file for general time output.";
	}


	return;
}

#endif /* _<GENERALHEADERH_ */
