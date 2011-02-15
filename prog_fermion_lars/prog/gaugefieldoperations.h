#ifndef _GAUGEFIELDOPSH_
#define _GAUGEFIELDOPSH_

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "operations.h"
#include "input.h"
#include "readgauge.h"
#include "random.h"
#include "use_timer.h"


void print_info_source(sourcefileparameters* params);

hmc_error init_gaugefield(hmc_gaugefield* gaugefield, inputparameters* parameters, usetimer* timer);

#endif
