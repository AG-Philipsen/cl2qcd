/** @file
 *
 * Everything required by heatbath's main()
 */
#ifndef _TKKAPPAH_
#define _TKKAPPAH_
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
#include "gaugefield_k.h"

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

void print_info(inputparameters* params, ostream* os)
{
	*os << "## **********************************************************\n";
	*os << "## Compile time parameters:\n";
	*os << "## NSPACE:  " << NSPACE << '\n';
	*os << "## NTIME:   " << NTIME << '\n';
	*os << "## NDIM:    " << NDIM << '\n';
	*os << "## NCOLOR:  " << NC << '\n';
	*os << "## NSPIN:   " << NSPIN << '\n';
	*os << "##" << '\n';
	*os << "## Run time parameters:\n";
	*os << "## beta  = " << params->get_beta() << '\n';
	*os << "## prec  = " << params->get_prec() << '\n';
	*os << "## thermsteps     = " << params->get_thermalizationsteps() << '\n';
	*os << "## heatbathsteps  = " << params->get_heatbathsteps() << '\n';
	*os << "## overrelaxsteps = " << params->get_overrelaxsteps() << '\n';
	*os << "##" << '\n';
	if (params->get_startcondition() == START_FROM_SOURCE) {
		*os << "## sourcefile = ";
		string sf=params->sourcefile;
		*os << sf << '\n';
	}
	if (params->get_startcondition() == COLD_START) {
		*os << "## cold start\n";
	}
	if (params->get_startcondition() == HOT_START) {
		*os << "## hot start\n";
	}
	*os << "## **********************************************************\n";
	*os << std::endl;
	return;
}

#endif /* _TKKAPPAH_ */
