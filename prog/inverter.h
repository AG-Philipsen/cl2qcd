/** @file
 *
 * Everything required by inverter's main()
 */
#ifndef _INVERTERH_
#define _INVERTERH_

//this includes header used by all executables
#include "general_header.h"

//should only be included in main prog
#include "gaugefield_inverter.h"
#include "gaugefield_hybrid.h"
string const version = "0.1b";

usetimer solver_timer;

void print_solver_profiling(std::string filename)
{
	uint64_t avg_time = 0.;
	int calls_total = solver_timer.getNumMeas();
	uint64_t time_total = solver_timer.getTime();
	//check if kernel has been called at all
	if(calls_total != 0 && time_total != 0) {
		avg_time = (uint64_t) ( ( (float) time_total ) / ((float) calls_total) );
	}
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) File_Exception(filename.c_str());
	out.width(32);
	out.precision(15);
	out << "## solver time [mus]:\ttime\tcalls\tavg" << endl;
	out << "\t" << time_total << "\t" << calls_total << "\t" << avg_time << std::endl;
	out.close();
	
	//logger.info() << "## solver time [mus]:\ttime\tcalls\tavg" ;
	//logger.info() << "\t" << time_total << "\t" << calls_total << "\t" << avg_time;
	return;
}


#endif /* _INVERTERH_ */

