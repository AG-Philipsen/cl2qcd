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
#include <iostream>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_geometry.h"
#include "host_input.h"
#include "host_readgauge.h"
#include "host_random.h"
#include "host_update_heatbath.h"
#include "host_use_timer.h"
#include "gaugefield.h"
#include "opencl.h"
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

//couple of timers
usetimer totaltime;
usetimer inittime;
usetimer polytime;
usetimer plaqtime;
usetimer updatetime;
usetimer overrelaxtime;
usetimer copytime;

//to save gaugeobservables
hmc_float plaq, splaq, tplaq;
hmc_complex pol;

void print_hello(char* name)
{
	std::cout << "This is heatbath program, " << name << endl;
	return;
}

void print_info(inputparameters* params)
{
	std::cout << "**********************************************************\n";
	std::cout << "Compile time parameters:\n";
	std::cout << "NSPACE:  " << NSPACE << '\n';
	std::cout << "NTIME:   " << NTIME << '\n';
	std::cout << "NDIM:    " << NDIM << '\n';
	std::cout << "NCOLOR:  " << NC << '\n';
	std::cout << "NSPIN:   " << NSPIN << '\n';
	std::cout << '\n';
	std::cout << "Run time parameters:\n";
	std::cout << "beta  = " << params->get_beta() << '\n';
	std::cout << "prec  = " << params->get_prec() << '\n';
	std::cout << "thermsteps     = " << params->get_thermalizationsteps() << '\n';
	std::cout << "heatbathsteps  = " << params->get_heatbathsteps() << '\n';
	std::cout << "overrelaxsteps = " << params->get_overrelaxsteps() << '\n';
	std::cout << '\n';
	if (params->get_startcondition() == START_FROM_SOURCE) {
		std::cout << "sourcefile = ";
		params->display_sourcefile();
		std::cout << '\n';
	}
	if (params->get_startcondition() == COLD_START) {
		std::cout << "cold start\n";
	}
	if (params->get_startcondition() == HOT_START) {
		std::cout << "hot start\n";
	}
	std::cout << "**********************************************************\n";
	std::cout << std::endl;
	return;
}

#endif /* _HEATBATH_ */
