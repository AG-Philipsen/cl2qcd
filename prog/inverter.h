/** @file
 *
 * Everything required by inverter's main()
 */
#ifndef _INVERTERH_
#define _INVERTERH_
//should only be included in main prog

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
#include "host_random.h"
#include "host_use_timer.h"
#include "gaugefield.h"
#include "gaugefield_inversion.h"
#include "logger.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif

string const version = "0.1b";
using namespace std;

//global random number generator
Random rnd (seed);

//couple of timers
usetimer total_timer;
usetimer init_timer;
usetimer perform_timer;

usetimer copy_to_from_dev_timer;
usetimer copy_on_dev_timer;
usetimer solver_timer;

void inverter_time_output(usetimer * total, usetimer * init_timer, usetimer * perform_timer, usetimer * copy_to_from_dev_timer, usetimer* copy_on_dev_timer, usetimer * solver_timer){

	uint64_t totaltime = (*total).getTime();

	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time
	
	uint64_t init_time = (*init_timer).getTime();
	uint64_t perform_time = (*perform_timer).getTime();
	uint64_t copy1_time = (*copy_on_dev_timer).getTime();
	uint64_t copy2_time = (*copy_to_from_dev_timer).getTime();
	uint64_t solvertime = (*solver_timer).getTime();

	int copy1_steps;
	int copy2_steps;
	int solver_steps;
	
	uint64_t copy1_avgtime;
	uint64_t copy2_avgtime;
	uint64_t solver_avgtime;
	
	copy2_steps = (*copy_on_dev_timer).getNumMeas();
	copy2_steps = (*copy_to_from_dev_timer).getNumMeas();
	solver_steps = (*solver_timer).getNumMeas();
	
	copy2_avgtime = divide(copy2_time, copy2_steps);
	copy1_avgtime = divide(copy1_time, copy1_steps);
	solver_avgtime = divide(solvertime, solver_steps);

	logger.trace() << "*******************************************************************";
	logger.trace() << "Times [mus]:";
	logger.trace() << "*******************************************************************";
	logger.trace() << "Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
	logger.trace() << "Total:\t" << setfill(' ') << setw(12) << totaltime;
	logger.trace() << "Init.:\t" << setfill(' ') << setw(12) << init_time << '\t'<< fixed << setw(5) << setprecision(1) << percent(init_time, totaltime) ;
	logger.trace() << "Perf.:\t" << setfill(' ') << setw(12) << perform_time << '\t'<< fixed << setw(5) << setprecision(1) << percent(perform_time, totaltime) ;
	logger.trace() << "*******************************************************************";
	logger.trace() << "Other:\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg"<< '\t' << setw(5) << "perc";
	logger.trace() << "Solve:\t" << setfill(' ') << setw(12) << solvertime << '\t' << setw(12) << solver_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(solvertime, totaltime);
	logger.trace() << "Copy1:\t" << setfill(' ') << setw(12) << copy1_time << '\t' << setw(12) << copy1_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(copy1_time, totaltime);
	logger.trace() << "Copy2:\t" << setfill(' ') << setw(12) << copy2_time << '\t' << setw(12) << copy2_avgtime << '\t'<< fixed << setw(5) << setprecision(1) << percent(copy2_time, totaltime);
	logger.trace() << "*******************************************************************";

	logger.trace() << "No output to file implemented yet...";
	/** @todo output to file is not implemented */
// 	ofstream out;
// 	stringstream str_filename;
// 	//CP: this H is for heatbath benchmarking, it has to be replaced meaningfully for other occasions
// 	str_filename << "time_";
// #ifdef _USEGPU_
// 	str_filename << "G_";
// #else
// 	str_filename << "C_";
// #endif
// #ifdef _USEDOUBLEPREC_
// 	str_filename << "D_";
// #else
// 	str_filename << "S_";
// #endif
// #ifdef _RECONSTRUCT_TWELVE_
// 	str_filename << "R_";
// #else
// 	str_filename << "N_";
// #endif
// 
// 	str_filename << "inverter";
// 	out.open(str_filename.str().c_str(), fstream::app);
// 	if (out.is_open()) {
// 		//output:
// 		//NTIME   NSPACE   VOL4D   totaltime   inittime  solvertime copytime singletime   scalarproducttime   latime  Mtime Mdiagtime dslashtime(all times average time-measurement)
// 		out <<
// 		    NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t"
// 		    << init_ferm << "\t" << solvertime << "\t" << copy_ferm_avgtime << "\t" << single_ferm_avgtime << "\t" << scalprod_avgtime << "\t" << la_avgtime << "\t" << M_avgtime << "\t" << Mdiag_avgtime << "\t" << dslash_avgtime << endl;
// 		out.close();
// 
// 	} else logger.error() << "Unable to open file for output";
	return;
}

#endif /* _INVERTERH_ */
