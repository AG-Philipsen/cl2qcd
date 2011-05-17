/** @file
 * Observables for fermions
 */

#ifndef _FERMIONOBSERVABLESH_
#define _FERMIONOBSERVABLESH_

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_geometry.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_solver.h"

/**
 * Specify whether to use even-odd preconditioning.
 */
extern int const use_eo;

/**
 * Calculates charged pion correlator
 * @param[in] parameters provides parameters needed
 * @param[in] gaugefield input gaugefield
 * @todo CP: the output is at the moment just printed to screen. this has to change at some point!
 */
hmc_error simple_correlator(inputparameters * parameters, hmc_gaugefield* gaugefield);


#endif /* _FERMIONOBSERVABLESH_ */
