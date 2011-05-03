/** @file
 * Gauge field handline
 */

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
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_input.h"
#include "host_readgauge.h"
#include "host_use_timer.h"
#include "host_writegaugefield.h"
#include "host_gaugeobservables.h"

extern string const version;

void print_info_source(sourcefileparameters* params);

/**
 * Initialize the gaugefield as specified by the input file.
 *
 * @param[out]    gaugefield The gaugefield to be initialized
 * @param[in]     parameters The parsed input file
 * @param[in,out] timer      Timer for logging of execution time
 * @return Error code as defined in hmcerrs.h:
 *         @li HMC_XMLERROR     if file to read gaugefield from cannot be parsed
 *         @li HMC_INVALIDVALUE if no valid start condition has been selected
 *         @li HMC_SUCCESS otherwise
 */
hmc_error init_gaugefield(hmc_gaugefield* gaugefield, inputparameters* parameters, usetimer* timer);

hmc_error save_gaugefield(hmc_gaugefield* gaugefield, inputparameters* parameters, int number);

#endif /* _GAUGEFIELDOPSH_ */
