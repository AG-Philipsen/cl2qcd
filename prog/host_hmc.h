/** @file
 * Implementation of steps of the HMC
 */

#ifndef _HOST_HMCH_
#define _HOST_HMCH_

#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_random.h"
#include "host_operations_complex.h"
#include "host_operations_matrix.h"
#include "host_operations_su3matrix.h"
#include "host_operations_gaugemomentum.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
#include "host_update_heatbath.h"
#include "host_gaugeobservables.h"
#include "host_use_timer.h"
#include "host_solver.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <cmath>

/**
 * Calculate (Wilson) Gauge-Action /f$S_{Gauge} = \ldots/f$
 * @param[in] field input gaugefield
 * @param[in] beta parameter of action
 * @return Value of the action (hmc_float)
 * @todo needs testing
 * @todo implement improvement terms
 * @todo insert latex formula here
 */
hmc_float s_gauge(hmc_gaugefield * field, hmc_float beta); 

/**
 * Calculate Hamiltonian /f$ H = S_{Gauge} + S_{Fermion}  + 1/2*P^2 /f$
 * @param[in] field input gaugefield
 * @param[in] beta parameter for /f$S_{Gauge}/f$
 * @param[in] p input gauge momentum
 * @param[in] phi parameter for /f$S_{Fermion}/f$ (only if _FERMIONS_ is set)
 * @param[in] MdaggerMphi parameter for /f$S_{Fermion}/f$ (only if _FERMIONS_ is set)
 * @return Value of the hamiltonian (hmc_complex)
 * @todo check the return values. They are supposed to be real!! If that is the case one can change the return argument
 * @todo needs testing
 */
hmc_complex hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_gauge_momentum * p 
#ifdef _FERMIONS_
	,hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi
#endif
	); 

/**
 * Perform Metropolis-Step.
 * If it has to be, performs the change old->new.
 * 
 * @param[in] rndnumber Random Number used for the actual Metropolis-Step in the end
 * @param[in] beta parameter for Hamiltonian
 * @param[in] phi parameter for Hamiltonian (only if _FERMIONS_ is set)
 * @param[in] MdaggerMphi parameter for Hamiltonian (only if _FERMIONS_ is set)
 * @param[in] field Old Gauge Configuration
 * @param[in] p Old Gauge Momentum
 * @param[in] new_field New Gauge Configuration
 * @param[in] new_p New Gauge Momentum
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 * @todo CP: if it works, one should replace the return value with 0 or 1, depending on wether the new config was accepted or not!!
 */
hmc_error metropolis(hmc_float rndnumber, hmc_float beta
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi
#endif
	, hmc_gaugefield * field, hmc_gauge_momentum * p, hmc_gaugefield * new_field, hmc_gauge_momentum * new_p);

/**
 * Molecular Dynamics Update of the Gauge Momenta using the Leapfrog-scheme:
 * /f[
 * p_{out} = p_{in} - \epsilon \text{force}_{in}
 * /f]
 * 
 * @param[in] eps Leapfrog stepsize
 * @param[in,out] p_inout input/output gauge momentum
 * @param[in] force_in input force
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 * 
 */
hmc_error md_update_gauge_momenta(hmc_float eps, hmc_gauge_momentum * p_inout, hmc_gauge_momentum * force_in); 

/**
 * Molecular Dynamics Update of the Gaugefield using the Leapfrog-scheme:
 * /f[
 * U_{out} = \exp\left(i \epsilon p_{in}\right) u_{in}
 * /f]
 * 
 * @param[in] eps Leapfrog stepsize
 * @param[in] p_in input gauge momentum
 * @param[in,out] u_in input/output gaugefield
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 * 
 */
hmc_error md_update_gaugefield(hmc_float eps, hmc_gauge_momentum * p_in, hmc_gaugefield * u_inout);

#ifdef _FERMIONS_
/**
 * Calculate fermion force for molecular dynamics. The Force is actually a gauge momentum vector.
 * It is assumed that the spinorfield /f$\phi/f$ has already been inverted.
 * @param[in] parameters parameters needed for the fermionmatrix
 * @param[in] field input gaugefield (constant in this function)
 * @param[in] phi input spinorfield
 * @param[in] phi_inv (inverted) input spinorfield
 * @param[out] out output gauge momentum
 * @return Error code as defined in hmcerrs.h
 * @todo implement actual code
 * @todo needs testing
 *
 */
hmc_error fermion_force(inputparameters * parameters, hmc_gaugefield * field, hmc_spinor_field * phi, hmc_spinor_field * phi_inv, hmc_gauge_momentum * out); 

/**
 * Molecular Dynamics Update of the Spinorfield:
 * /f[
 * \phi = M \xi,
 * /f]
 * with M being the fermionmatrix and /f$ \xi /f$ a gaussian vector
 * @param[in] in input spinorfield (gaussian)
 * @param[out] out output spinorfield
 * @param[in] field input gaugefield
 * @param[in] parameters all parameters, the function extract those needed for the fermionmatrix
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 * @todo check if it is M or Mdagger
 * @todo extract needed parameters inside the functions!
 * 
 */
hmc_error md_update_spinorfield(hmc_spinor_field * in, hmc_spinor_field * out, hmc_gaugefield * field, inputparameters * parameters);

/**
 * Calculate Fermion-Action /f$S_{Fermion} = \phi ( M^{\dagger}M \phi ) /f$
 
 * @param[in] phi input spinorfield /f$\phi/f$
 * @param[in] MdaggerMphi input spinorfield /f$M^{\dagger}M \phi/f$
 * @return Value of the action
 * @todo needs testing
 * @todo check if return value is always real. If so, change the return argument.
 *
 */
hmc_complex s_fermion(hmc_spinor_field * phi, hmc_spinor_field * MdaggerMphi);
#endif

/**
 * Calculates the Gauge-Force. The Force is actually a gauge momentum vector.
 *
 * @param[in] parameters input parameters
 * @param[in] field input gaugefield
 * @param[out] out output gauge momentum
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 *
 */
hmc_error gauge_force(inputparameters * parameters, hmc_gaugefield * field, hmc_gauge_momentum * out);

/**
 * Calculates the force for the molecular dynamics.
 * It is assumed that phi has already been inverted
 *
 * @param[in] parameters input parameters 
 * @param[in] field input gaugefield
 * @param[in] phi input spinorfield (only if _FERMIONS_ is set)
 * @param[in] phi_inv (inverted) input spinorfield (only if _FERMIONS_ is set)
 * @param[out] out output gauge momentum
 * @return Error code as defined in hmcerrs.h
 * @todo needs testing
 */
hmc_error force(inputparameters * parameters, hmc_gaugefield * field
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
#endif
	, hmc_gauge_momentum * out);

/**
 * Performs the Leapfrog discretised Molecular Dynamics as in Gattringer/Lang, QCD on the Lattice, 8.2, p 197.
 * @param[in] parameters input parameters
 * @param[in] u_in input gaugefield
 * @param[in] p_in input gauge momentum
 * @param[in] phi input spinorfield (if _FERMIONS_ is on)
 * @param[out] u_out output gaugefield
 * @param[out] p_out output gauge momentum
 * @param[out] phi_inv inverted spinorfield (if _FERMIONS_ is on). This is used again outside of this function.
 * @return Error code as defined in hmcerrs.h
 *
 * @todo lateron, a multi-step alg. with different stepsizes for gauge and fermion force should be implemented
 * @todo see code itself for more points
 */
hmc_error leapfrog(inputparameters * parameters, hmc_gaugefield * u_in, hmc_gauge_momentum * p_in
#ifdef _FERMIONS_
	, hmc_spinor_field * phi, hmc_spinor_field * phi_inv
#endif
	, hmc_gaugefield * u_out, hmc_gauge_momentum * p_out
	); 


#endif /* _HOST_HMCH_ */
