/** @file
 * Implementation of steps of the HMC
 */

#ifndef _HMCH_
#define _HMCH_

#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_random.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
#include "host_geometry.h"
#include "host_gaugeobservables.h"
#include "host_use_timer.h"
#include "host_solver.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <cmath>

hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta); //TODO CP: not tested

hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p 
#ifdef _FERMIONS_
	,hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi
#endif
	); //TODO CP: not tested

hmc_error metropolis(hmc_float rndnumber, hmc_float beta
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi
#endif
	, hmc_gaugefield * field, hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p); // TODO SL: not tested

hmc_error md_update_gauge_momenta(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gauge_momentum * force_in, hmc_gauge_momentum * p_out); //TODO CP: not tested
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_in); //TODO CP: not tested

#ifdef _FERMIONS_
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_gauge_momentum * out); //TODO CP: not tested
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters); //TODO CP: not tested
hmc_error generate_gaussian_spinorfield(hmc_spinor_field * out); //TODO CP: not tested
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi); //TODO CP: not tested
#endif

hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out); //TODO CP: not tested
hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_gauge_momentum * out); //TODO CP: not tested

hmc_error force(inputparameters * parameters, hmc_gaugefield * field
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
#endif
	, hmc_gauge_momentum * out); //TODO CP: not tested

hmc_error leapfrog(inputparameters * parameters, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in
#ifdef _FERMIONS_
	, hmc_spinor_field * phi
#endif
	, hmc_gaugefield * u_out, hmc_gauge_momentum * p_out
#ifdef _FERMIONS_
	, hmc_spinor_field * phi_inv
#endif
	); //TODO CP: not tested

hmc_error construct_3x3_combination(hmc_float beta_0, hmc_float gamma_0, hmc_float beta[], hmc_float gamma[], hmc_3x3matrix out); //TODO CP: not tested
hmc_error build_su3matrix_by_exponentiation(hmc_algebraelement in, hmc_su3matrix out, hmc_float epsilon); //TODO CP: not tested

#endif /* _HMCH_ */
