#ifndef _HMCH_
#define _HMCH_
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_random.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
#include "host_gaugeobservables.h"
#include "host_use_timer.h"
#include "host_solver.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <cmath>

hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta); //CP: not tested
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi); //CP: not tested
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi); //CP: not tested
hmc_error metropolis(hmc_float rndnumber, hmc_float beta, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi, hmc_gaugefield * field,
					 hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p); // SL: not tested

#endif 