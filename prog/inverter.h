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
#include "gaugefield.h"
#include "gaugefield_inversion.h"

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
	std::cout << "This is the fermion solver, " << name << endl;
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

	if(params->get_fermact()==WILSON) {
	  *os<<  "## fermion action: unimproved Wilson"<<'\n';
	  *os << "## kappa  = "<<params->get_kappa()<< '\n';
	}
	if(params->get_fermact()==TWISTEDMASS) {
	  *os<<  "## fermion action: twisted mass Wilson"<<'\n';
	  *os << "## kappa  = "<<params->get_kappa()<< '\n';
	  *os << "## mu     = "<<params->get_mu()<< '\n';
	}
	if(params->get_fermact()==CLOVER) {
	  *os<<  "## fermion action: clover Wilson"<<'\n';
	  *os << "## kappa  = "<<params->get_kappa()<< '\n';
	  *os << "## csw    = "<<params->get_csw()<< '\n';
	}
	*os << "## prec  = " << params->get_prec() << '\n';
	*os << "##" << '\n';

	if (params->get_startcondition() == START_FROM_SOURCE) {
		*os << "## sourcefile = ";
		string sf = params->sourcefile;
		*os << sf << '\n';
	}
	if (params->get_startcondition() == COLD_START) {
		*os << "## WARNING: cold start - no configuration read\n";
	}
	if (params->get_startcondition() == HOT_START) {
		*os << "## WARNING: hot start - no configuration read\n";
	}
	*os << "## **********************************************************\n";
	*os << std::endl;
	return;
}

#endif /* _INVERTERH_ */
