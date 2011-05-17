/** @file
 *
 * Everything required by main()
 */

#ifndef _HMCH_
#define _HMCH_
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
#include "host_operations_matrix.h"
#include "host_operations_su3matrix.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_geometry.h"
#include "host_testing.h"
#include "host_gaugeobservables.h"
#include "host_gaugefieldoperations.h"
#include "host_input.h"
#include "host_readgauge.h"
#include "host_random.h"
#include "host_update_heatbath.h"
#include "host_use_timer.h"
#include "host_writegaugefield.h"
#include "opencl.h"
#include "host_hmc.h"
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

string const version = "0.1";

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
#ifdef _FERMIONS_
usetimer inittimer;
usetimer singletimer;
usetimer Mtimer;
usetimer dslashtimer;
usetimer Mdiagtimer;
usetimer copytimer;
usetimer scalarprodtimer;
usetimer latimer;
usetimer solvertimer;
#endif
#ifdef _USEHMC_
usetimer hmctimer;
usetimer leapfrogtimer;
usetimer hmcinittimer;
usetimer metropolistimer;
#endif


//to save gaugeobservables
hmc_float plaq, splaq, tplaq;
hmc_complex pol;

#ifdef _PERFORM_BENCHMARKS_
char * benchmark_id;
#endif

void print_hello(char* name)
{
	std::cout<<"This is hmc program "<<name<<", version "<<version<<"."<<endl;
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
	std::cout << "kappa = " << params->get_kappa() << '\n';
	std::cout << "mu    = " << params->get_mu() << '\n';
	std::cout << "csw   = " << params->get_csw() << '\n';
	std::cout << "beta  = " << params->get_beta() << '\n';
	std::cout << "CGmax = " << params->get_cgmax() << '\n';
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

#endif /* _HMCH_ */
