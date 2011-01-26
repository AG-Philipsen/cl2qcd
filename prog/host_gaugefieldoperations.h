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
#include "host_operations.h"
#include "host_input.h"
#include "host_readgauge.h"
#include "host_use_timer.h"


void print_info_source(sourcefileparameters* params);

hmc_error init_gaugefield(hmc_gaugefield* gaugefield, inputparameters* parameters, usetimer* timer);

#endif
