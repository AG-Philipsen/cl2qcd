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

// SL: the following flags on how to evaluate group = exp(i*epsilon*algebra) should maybe be moved elsewhere?
#define EXPONENTIATE_ALGEBRA_ORDER_2
	// #define EXPONENTIATE_ALGEBRA_ORDER_3
	// #define EXPONENTIATE_ALGEBRA_ALL_ORDERS
// Leave uncommented exactly one of the three above!

// SL: also these utilitied shall me placed somewhere else?
// Definition of numeric constants for the symmetric structure constants d_ijk of su(3)
#define F_1_2   (static_cast<hmc_float>(0.5))					// 1/2
#define F_1_2S3 (static_cast<hmc_float>(0.288675134594813))		// 1/(2*sqrt(3))
#define F_1_S3  (static_cast<hmc_float>(0.577350269189626))		// 1/sqrt(3)

hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta); //CP: not tested
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi); //CP: not tested
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi); //CP: not tested
hmc_error metropolis(hmc_float rndnumber, hmc_float beta, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi, hmc_gaugefield * field,
					 hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p); // SL: not tested
hmc_error md_update_gauge_momenta(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gauge_momentum * force_in, hmc_gauge_momentum * p_out);
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_in);
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters);
hmc_error generate_gaussian_spinorfield(hmc_spinor_field * out);
hmc_error generate_gaussian_gauge_momenta(hmc_gauge_momentum * out);
hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_gauge_momentum * out);
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * phi_inv,
					hmc_gauge_momentum * out);
hmc_error force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_gauge_momentum * out,
					hmc_spinor_field * phi_inv);
hmc_error leapfrog(inputparameters * parameters, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in, hmc_spinor_field * phi,
					hmc_gaugefield * u_out, hmc_gauge_momentum * p_out, hmc_spinor_field * phi_inv);
hmc_error construct_3x3_combination(hmc_float beta_0, hmc_float gamma_0, hmc_float beta[], hmc_float gamma[], hmc_3x3matrix out);
hmc_error build_su3matrix_by_exponentiation(hmc_algebraelement in, hmc_su3matrix out, hmc_float epsilon);
#endif 


