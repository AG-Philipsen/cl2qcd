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

hmc_error simple_correlator(hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax);

#endif