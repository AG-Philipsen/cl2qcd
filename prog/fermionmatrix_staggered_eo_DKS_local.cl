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
  *       get_staggered_phase \internal(see operations_staggered.cl).\endinternal
  * \par
  * @note If we chose to impose boundary conditions modifying staggered phases at the end of the
  *       lattice in each direction, then we would make each staggered phase appear EXACTELY
  *       next to the link, and if the link is dagger, we would take the complex coniugate
  *       of the staggered phase (that would be complex in general due to the modification).
  * 
  * @todo If a chemical potential is introduced, this kernel has to be modified!
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
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	if(dir == TDIR){
	  bc_tmp.re = TEMPORAL_RE;
	  bc_tmp.im = TEMPORAL_IM;
	}else{
	  bc_tmp.re = SPATIAL_RE;
	  bc_tmp.im = SPATIAL_IM;
	}
	//this is used to take into account the staggered phase
	hmc_float eta;
	out_tmp = set_su3vec_zero();
	
	//go through the different directions
	///////////////////////////////////
	// mu = +dir
	///////////////////////////////////
	idx_neigh = get_neighbor_from_st_idx(idx_arg, dir);
	//transform normal indices to eoprec index
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);
	plus = get_su3vec_from_field_eo(in, nn_eo);
	U = getSU3(field, get_link_idx(dir, idx_arg));
	//Thanks to the variables passed to this kernel I can write all directions in few lines:
	//chi=U*plus
	chi=su3matrix_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_acc(out_tmp,chi);
	
	//The 3 lines above are equivalent to the following code
	/*
	if(dir == XDIR) {
		//chi=U*plus
		chi=su3matrix_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_acc(out_tmp,chi);
		
	} else if(dir == YDIR) {
		//chi=U*plus
		chi=su3matrix_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_acc(out_tmp,chi);
	  
	} else if(dir == ZDIR) {
		//chi=U*plus
		chi=su3matrix_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_acc(out_tmp,chi);
	  
	} else { // dir == TDIR
		//chi=U*plus
		chi=su3matrix_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_acc(out_tmp,chi);

	}
	*/
	
	///////////////////////////////////
	// mu = -dir
	///////////////////////////////////
	idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir);
	//transform normal indices to eoprec index
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);
	plus = get_su3vec_from_field_eo(in, nn_eo);
	U = getSU3(field, get_link_idx(dir, idx_neigh));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	if(dir == TDIR){
	  bc_tmp.re = TEMPORAL_RE;
	  bc_tmp.im = MTEMPORAL_IM;
	}else{
	  bc_tmp.re = SPATIAL_RE;
	  bc_tmp.im = MSPATIAL_IM;
	}
	//Thanks to the variables passed to this kernel I can write all directions in few lines:
	//chi=U^dagger*plus
	chi=su3matrix_dagger_times_su3vec(U,plus);
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_dim(out_tmp,chi);
	//multiply by the factor 1/2*eta_dir that appears at the beginning of D_KS
	//if dir==XDIR then eta_x is 1 and the multiplication is unnecessary
	if(dir == XDIR)
		out_tmp = su3vec_times_real(out_tmp, F_1_2); //AS: I'm not sure that here F_1_2 is an hmc_float
	else{
		eta=0.5*get_staggered_phase(idx_arg.space,idx_arg.time,dir);
		out_tmp = su3vec_times_real(out_tmp, eta);
	}
	//The 9 lines of code above are equivalent to the following code
	/*
	if(dir == XDIR) {
		//chi=U^dagger*plus
		chi=su3matrix_dagger_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_dim(out_tmp,chi);
	} else if(dir == YDIR) {
		//chi=U^dagger*plus
		chi=su3matrix_dagger_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_dim(out_tmp,chi);
		//multiply by the factor 1/2*eta_y that appears at the beginning of D_KS
		eta=0.5*get_staggered_phase(idx_arg.space,idx_arg.time,dir);
		out_tmp = su3vec_times_real(out_tmp, eta); 
		
	} else if(dir == ZDIR) {
		//chi=U^dagger*plus
		chi=su3matrix_dagger_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_dim(out_tmp,chi);
		//multiply by the factor 1/2*eta_z that appears at the beginning of D_KS
		eta=0.5*get_staggered_phase(idx_arg.space,idx_arg.time,dir);
		out_tmp = su3vec_times_real(out_tmp, eta); 
		  
	} else { // TDIR
		//chi=U^dagger*plus
		chi=su3matrix_dagger_times_su3vec(U,plus);
		chi=su3vec_times_complex(chi,bc_tmp);
		out_tmp=su3vec_dim(out_tmp,chi);
		//multiply by the factor 1/2*eta_t that appears at the beginning of D_KS
		eta=0.5*get_staggered_phase(idx_arg.space,idx_arg.time,dir);
		out_tmp = su3vec_times_real(out_tmp, eta); 
	}
	*/

	return out_tmp;
}

