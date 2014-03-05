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
* Staggered (local) D_KS operator that is called from within kernels \b with EVEN-ODD preconditioning
* \internal (for example, see fermionmatrix_staggered.cl)
*/

/**
  * This kernel is nothing but the local D_KS working on a particular link in a specific direction, 
  * taking into account even-odd preconditioning. Indeed, here there are no conceptual aspects to
  * be implemented to take into account eo-prec. The difference between this kernel and D_KS_local
  * is that here the coordinates of the neighbors have to be transformed into an eoprec index and then
  * the functions get_su3vec_from_field_eo must be used.
  * \internal (see spinorfield_staggered_eo.cl for these 2 functions) \endinternal  
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
  *  @param in The input staggered field on half lattice (either on even or odd sites)
  *  @param field The links configuration
  *  @param idx_arg The superindex of the site where the output field is returned
  *  @param dir The direction in which D_KS works
  * @return The field in the idx_arg site is returned. Remark that if the "in" staggered
  *         field is on even sites then idx_arg will be an odd site and vice-versa.
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
  * \par 
  * @note If an imaginary chemical potential is used, then the links in time direction have to
  *       be multiplied by the phases exp(i\mu) and exp(-i\mu). Actually we multiply temporal
  *       links only by exp(i\mu) because then, backward in time, we use U dagger.
  */
su3vec D_KS_eo_local(__global const staggeredStorageType * const restrict in, __global const Matrixsu3StorageType * const restrict field, const st_idx idx_arg, const dir_idx dir)
{
	//this is used to save the idx of the neighbors
	st_idx idx_neigh;
	//this is used to transform the idx of the neighbors, that is a superindex (between 0 to VOL4D-1),
	//to a superindex of type even-odd (between 0 and VOL4D/2-1)
	site_idx nn_eo;
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
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);//transform normal indices to eoprec index
	plus = get_su3vec_from_field_eo(in, nn_eo);
	U = getSU3(field, get_link_idx(dir, idx_arg));
#ifdef _CP_IMAG_
	//Simplest code, if low performance try something like (dir==TDIR)*cpi_tmp to avoid the if.
	//Actually one could also think to include the imaginary chemical potential
	//in the staggered phases as done for the boundary conditions. In this case one
	//should move cpi_tmp to the file operations_staggered.cl
	if(dir == TDIR){ 
	  hmc_complex cpi_tmp = {COSCPI, SINCPI};
	  U = multiply_matrixsu3_by_complex(U, cpi_tmp);
	}
#endif
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
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);//transform normal indices to eoprec index
	plus = get_su3vec_from_field_eo(in, nn_eo);
	U = getSU3(field, get_link_idx(dir, idx_neigh));
#ifdef _CP_IMAG_
	//Simplest code, if low performance try something like (dir==TDIR)*cpi_tmp to avoid the if.
	//Actually one could also think to include the imaginary chemical potential
	//in the staggered phases as done for the boundary conditions. In this case one
	//should move cpi_tmp to the file operations_staggered.cl
	if(dir == TDIR){
	  hmc_complex cpi_tmp = {COSCPI, SINCPI};
	  U = multiply_matrixsu3_by_complex(U, cpi_tmp);
	}
#endif
	//chi=U^dagger * plus
	chi = su3matrix_dagger_times_su3vec(U, plus);
	eta_mod = get_modified_stagg_phase(idx_arg.space, dir);
	eta_mod.re *= 0.5; //the factors 0.5 is to take into 
	eta_mod.im *= 0.5; //account the factor in front of D_KS
	chi = su3vec_times_complex_conj(chi, eta_mod); //here conj is crucial for BC that are next to a U^dagger
	
	out_tmp=su3vec_dim(out_tmp,chi);
	
	return out_tmp;
}
