/** @file
 *
 * Everything required by inverter's main()
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
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "gaugefield.h"
#include "gaugefield_inversion.h"
#include "gaugefield_hmc.h"
#include "logger.hpp"
#include "types_hmc.h"

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

usetimer ferm_inittime;

usetimer copytimer;
usetimer singletimer;
usetimer Mtimer;
usetimer scalarprodtimer;
usetimer latimer;
usetimer dslashtimer;
usetimer Mdiagtimer;
usetimer solvertimer;


//to save gaugeobservables
hmc_float plaq, splaq, tplaq;
hmc_complex pol;


#endif /* _HMCH_ */
