/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
* Staggered (local) D_KS operator that is called from within kernels
* \internal (for example, see fermionmatrix_staggered_M.cl)
*/

/**
  * This kernel does exactely what is done in the four kernels dslash_local_i (see below in this file).
  * It is nothing but the local D_KS working on a particular link in a specific direction.
  * The expression of D_KS for a specific (couple of) site(s) and a specific direction is
  *  \f[
     (D_{KS})_{n,m,\mu}=\frac{1}{2} \eta_\mu(n)\Bigl[U_\mu(n)\,\delta_{n+\hat\mu,m} - U^\dag_\mu(n-\hat\mu)\,\delta_{n-\hat\mu,m}\Bigr]
     \f]
  * This function returns the value of the field (D_KS*in) at the site idx_arg [so in the
  * function only values of the field "in" in the nextneighbour of idx_arg will be needed]:
  *  \f[
     \bigl[(D_{KS})_\mu\cdot \text{\texttt{in}}\bigr]_n=\frac{1}{2}\eta_\mu(n) \Bigl[U_\mu(n) \cdot\text{\texttt{in}}_{n+\hat\mu} - U^\dag_\mu(n-\hat\mu)\cdot\text{\texttt{in}}_{n-\hat\mu}\Bigr]
    \f]
  *
  * The variables passed to the kernel are:
  *  @param in The input staggered field on the whole lattice
  *  @param field The links configuration
  *  @param idx_arg The superindex of the site where the output field is returned
  *  @param dir The direction in which D_KS works
  * 
  * @note The staggered phases are included in this function with the help of the function
  *       get_modified_stagg_phase @internal(see operations_staggered.cl)@endinternal.
  * \par
  * @note Since we chose to impose boundary conditions as in the Wilson Code, we have to multiply
  *       each link by a proper phase (and each link dagger by the conjugated phase).
  *       This is done via including this BC phase in the staggered phase calculation making
  *       it a complex number. Then we have to multiply each link separately by the staggered
  *       phase, in order to take correctly the complex conjugated. This makes the code more
  *       symmetric and raises the performance.
  * 
  * @todo If a chemical potential is introduced, this kernel has to be modified!
  */
su3vec D_KS_local(__global const su3vec * const restrict in, __global const Matrixsu3StorageType * const restrict field, const st_idx idx_arg, const dir_idx dir)
{
	//this is used to save the idx of the neighbors
	st_idx idx_neigh;
	//this are used for the calculation
	su3vec out_tmp, plus, chi;
	Matrixsu3 U;

	//this is used to take into account the staggered phase and the BC-conditions
	hmc_complex eta_mod;
	
	out_tmp = set_su3vec_zero();
	
	//go through the different directions
	///////////////////////////////////
	// mu = +dir
	///////////////////////////////////
	idx_neigh = get_neighbor_from_st_idx(idx_arg, dir);
	plus = get_su3vec_from_field(in, idx_neigh.space, idx_neigh.time);
	U = getSU3(field, get_link_idx(dir, idx_arg));
	//chi=U*plus
	chi = su3matrix_times_su3vec(U, plus);
	eta_mod = get_modified_stagg_phase(idx_arg.space, dir);
	eta_mod.re *= 0.5; //the factors 0.5 is to take into 
	eta_mod.im *= 0.5; //account the factor in front of D_KS
	chi = su3vec_times_complex(chi, eta_mod);
	
	out_tmp = su3vec_acc(out_tmp, chi);
	
	///////////////////////////////////
	// mu = -dir
	///////////////////////////////////
	idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir);
	plus = get_su3vec_from_field(in, idx_neigh.space, idx_neigh.time);
	U = getSU3(field, get_link_idx(dir, idx_neigh));
	//chi=U^dagger * plus
	chi = su3matrix_dagger_times_su3vec(U, plus);
	eta_mod = get_modified_stagg_phase(idx_arg.space, dir);
	eta_mod.re *= 0.5; //the factors 0.5 is to take into 
	eta_mod.im *= 0.5; //account the factor in front of D_KS
	chi = su3vec_times_complex_conj(chi, eta_mod); //here conj is crucial for BC that are next to a U^dagger
	
	out_tmp = su3vec_dim(out_tmp, chi);
	
	return out_tmp;
}






#if 0
// At the moment we use the kernel above. Uncomment the following part of code if needed for some reason.

/** \e Local D_KS working on a particular link (n,t) in a specific direction. For the expression of D_KS
 *  for a specific (couple of) site(s) and a specific direction see D_KS_local kernel documentation.
 *  This function returns the value of the field (D_KS*in) at the site (n,t) [so in the
 *  function only values of the field "in" in the nextneighbour of (n,t) will be needed].
 * 
 * @note Again, the staggered phases are included in this function with the help of the function
 *       get_staggered_phase.
 * 
 * @todo Here there are some #ifdef about the chemical potential: they have been copied from the
 *       Wilson code. Therefore they MUST be checked when a chemical potential will be introduced.
 */

su3vec dslash_local_0(__global const su3vec * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t)
{
	su3vec out_tmp, plus, chi;
	int dir, nn;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_su3vec_zero();
	//this is used to take into account the staggered phase
	hmc_float eta;
	
	//go through the different directions
	///////////////////////////////////
	// mu = 0
	///////////////////////////////////
	dir = 0;
	///////////////////////////////////
	//mu = +0
	/////////////
	nn = get_neighbor_temporal(t);
	plus = get_su3vec_from_field(in, n, nn);
	U = getSU3(field, get_link_pos(dir, n, t));
	//if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
	bc_tmp.re = TEMPORAL_RE;
	bc_tmp.im = TEMPORAL_IM;
	///////////////////////////////////
	//chi=U*plus
	////////////////
	chi=su3matrix_times_su3vec(U,plus);
	//chi=su3vec_times_real(chi,0.5*get_staggered_phase(n,t,dir));
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_acc(out_tmp,chi);
	
	/////////////////////////////////////
	//mu = -0
	//////////////
	nn = get_lower_neighbor_temporal(t);
	plus = get_su3vec_from_field(in, n, nn);
	U = getSU3(field, get_link_pos(dir, n, nn));
	//if chemical potential is activated, U has to be multiplied by appropiate factor
	//this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
	//  (exp(iq)U)^+ = exp(-iq)U^+
	//as it should be
	//in the real case, one has to take exp(q) -> exp(-q)
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, MEXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp2 = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp2 );
#endif
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = TEMPORAL_RE;
	bc_tmp.im = MTEMPORAL_IM;
	/////////////////////////////////////
	//chi=U^dagger*plus
	////////////////////////
	chi=su3matrix_dagger_times_su3vec(U,plus);
	//chi=su3vec_times_real(chi,0.5*get_staggered_phase(n,nn,dir));
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_dim(out_tmp,chi);
	
	///////////////////////////////////
	//multiply by the factor 1/2*eta_t that appears at the beginning of D_KS
	///////////////////////////////////
	eta=0.5*get_staggered_phase(n,t,dir);
	out_tmp = su3vec_times_real(out_tmp, eta); 

	return out_tmp;
}

su3vec dslash_local_1(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t)
{
	su3vec out_tmp, plus, chi;
	int dir, nn;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_su3vec_zero();
	
	//CP: all actions correspond to the mu = 0 ones
	///////////////////////////////////
	// mu = 1
	///////////////////////////////////
	dir = 1;
	///////////////////////////////////
	// mu = +1
	//////////////
	nn = get_neighbor_spatial(n, dir);
	plus = get_su3vec_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, n, t));
	bc_tmp.re = SPATIAL_RE;
	bc_tmp.im = SPATIAL_IM;
	/////////////////////////////////
	//chi=U*plus
	//////////////
	chi=su3matrix_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_acc(out_tmp,chi);
	
	///////////////////////////////////
	//mu = -1
	/////////////
	nn = get_lower_neighbor_spatial(n, dir);
	plus = get_su3vec_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = SPATIAL_RE;
	bc_tmp.im = MSPATIAL_IM;
	///////////////////////////////////
	//chi=U^dagger*plus
	//////////////////////
	chi=su3matrix_dagger_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_dim(out_tmp,chi);
	
	///////////////////////////////////
	//multiply by the factor 1/2 that appears at the beginning of D_KS
	//Observe that in the x direction ALL staggered phases are +1 => no need to multiply
	///////////////////////////////////
	out_tmp = su3vec_times_real(out_tmp, F_1_2); //AS: I'm not sure that here F_1_2 is an hmc_float
		
	return out_tmp;
}

su3vec dslash_local_2(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t)
{
	su3vec out_tmp, plus, chi;
	int dir, nn;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_su3vec_zero();
	//this is used to take into account the staggered phase
	hmc_float eta;

	///////////////////////////////////
	//mu = 2
	///////////////////////////////////
	dir = 2;
	///////////////////////////////////
	//mu = +2
	////////////
	nn = get_neighbor_spatial(n, dir);
	plus = get_su3vec_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, n, t));
	bc_tmp.re = SPATIAL_RE;
	bc_tmp.im = SPATIAL_IM;
	///////////////////////////////////
	//chi=U*plus
	///////////////
	chi=su3matrix_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_acc(out_tmp,chi);

	///////////////////////////////////
	//mu = -2
	//////////////
	nn = get_lower_neighbor_spatial(n, dir);
	plus = get_su3vec_from_field(in,  nn, t);
	U = getSU3(field, get_link_pos(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = SPATIAL_RE;
	bc_tmp.im = MSPATIAL_IM;
	///////////////////////////////////
	//chi=U^dagger*plus
	////////////////////
	chi=su3matrix_dagger_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_dim(out_tmp,chi);
	
	///////////////////////////////////
	//multiply by the factor 1/2*eta_y that appears at the beginning of D_KS
	///////////////////////////////////
	eta=0.5*get_staggered_phase(n,t,dir);
	out_tmp = su3vec_times_real(out_tmp, eta);
	
	return out_tmp;
}

su3vec dslash_local_3(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t)
{
	su3vec out_tmp, plus, chi;
	int dir, nn;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_su3vec_zero();
	//this is used to take into account the staggered phase
	hmc_float eta;

	///////////////////////////////////
	//mu = 3
	///////////////////////////////////
	dir = 3;
	///////////////////////////////////
	//mu = +3
	///////////////
	nn = get_neighbor_spatial(n, dir);
	plus = get_su3vec_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, n, t));
	bc_tmp.re = SPATIAL_RE;
	bc_tmp.im = SPATIAL_IM;
	///////////////////////////////////
	//chi=U*plus
	///////////////////////////////////
	chi=su3matrix_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_acc(out_tmp,chi);

	///////////////////////////////////
	//mu = -3
	//////////////
	nn = get_lower_neighbor_spatial(n, dir);
	plus = get_su3vec_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = SPATIAL_RE;
	bc_tmp.im = MSPATIAL_IM;
	///////////////////////////////////
	//chi=U^dagger*plus
	///////////////////////////////////
	chi=su3matrix_dagger_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_dim(out_tmp,chi);
	
	///////////////////////////////////
	//multiply by the factor 1/2*eta_z that appears at the beginning of D_KS
	///////////////////////////////////
	eta=0.5*get_staggered_phase(n,t,dir);
	out_tmp = su3vec_times_real(out_tmp, eta); 
	
	return out_tmp;
}

#endif