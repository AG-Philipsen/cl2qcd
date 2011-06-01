/** @file
 * Operations on Gauge momentum.
 *
 * @todo A lot of these operations use pointers and in-place operation
 *       without need. These might prevent compiler optimizations.
 */
#ifndef _OPERATIONS_GAUGEMOMENTUMH_
#define _OPERATIONS_GAUGEMOMENTUMH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_random.h"
#include <cmath>

/**
 * Create a copy of the gauge momenta.
 *
 * @param[in] source The gauge momenta to copy
 * @param[out] dest The storage location for the copy
 * @return Error code as defined in hmcerrs.h
 */
hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest);
/**
 * Calculate the squarenorm of the gauge momenta.
 *
 * @param[in] in The gauge momenta to use.
 * @param[out] result The square norm
 * @return Error code as defined in hmcerrs.h
 */
hmc_error gaugemomenta_squarenorm(hmc_gauge_momentum * in, hmc_float * result);
/**
 * Set gaugemomenta to zero.
 *
 * @param[out] in The gauge momenta to set to zero.
 * @return Error code as defined in hmcerrs.h
 */
hmc_error set_zero_gaugemomenta(hmc_algebraelement2 * in);
//deprecated:
//hmc_error set_zero_gaugemomenta(hmc_gauge_momentum * in);


/**
 * Generates a gaussian distributed complex vector of length GAUGEMOMENTASIZE and variance 1.
 * @param[out] out output gauge momentum
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 *
 */
hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out); 

#endif
