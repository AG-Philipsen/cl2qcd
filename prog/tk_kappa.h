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
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_random.h"
#include "host_use_timer.h"
#include "gaugefield.h"
#include "gaugefield_heatbath.h"
#include "gaugefield_k.h"
#include "logger.hpp"

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

#endif /* _TKKAPPAH_ */
