/** @file
* Staggered fermionmatrix functions that are called from within kernels
* \internal (for exemple, see fermionmatrix_m_staggered.cl)
*/

/**
 * This functions returns the value of the field "in" at the site (n,t)
 *
 * @todo This function should go into an own file...
 */
su3vec get_su3vec_from_field(__global const su3vec * const restrict in, const int n, const int t)
{
	int pos = get_global_pos(n, t);
	su3vec out;
	out = in[pos];
	return out;
}

/**
 * This functions uploads the value of the field "out" at the site (n,t) with the value "in"
 *
 * @todo This function should go into an own file...
 */
void put_su3vec_to_field(const su3vec in, __global su3vec * const restrict out, const int n, const int t)
{
	int pos = get_global_pos(n, t);
	out[pos] = in;
}

/**
 * This function calculates the staggered phase lying at the site with spatial super-index
 * n and in direction dir.
 *
 * @note The whole site superindex is not needed because to calculate the staggered phase one
 *       never needs the time coordinate:
 *       \f[
             \left\{
             \begin{aligned}
                \eta_1(n)&=1 \\
                \eta_\mu(n)&=(-1)^{\sum_{\nu<\mu}n_\nu} \quad\mbox{if}\quad \mu\neq1
             \end{aligned}
             \right.
 *       \f]
 *       where is CRUCIAL to emphasize that
 *       \f[
             \begin{aligned}
               \mu=1 &\quad\mbox{means}\quad x=n_1 \\
               \mu=2 &\quad\mbox{means}\quad y=n_2 \\
               \mu=3 &\quad\mbox{means}\quad z=n_3 \\
               \mu=4 &\quad\mbox{means}\quad t=n_4 
             \end{aligned}
 *       \f]
 *
 * \internal
 * 
 * Strictly speaking we have:
 *
 *  In x direction --> eta_x=1
 *  In y direction --> eta_y=(-1)^(x)
 *  In z direction --> eta_z=(-1)^(x+y)
 *  In t direction --> eta_t=(-1)^(x+y+z)
 *
 * \endinternal
 *
 * @param n this is the spacial superindex of the site (see st_idx.space in operation_geometry.cl)
 * @param dir this is the direction of the staggered phase. To be automatically coherent with
 *             the choice of labels made in operation_geometry.cl, it can be YDIR, ZDIR or TDIR.
 *      
 * @todo This function must be moved to a file containing all staggered utilities.
 */
 
int get_staggered_phase(const int n, const int dir)
{
	coord_spatial coord;
	coord=get_coord_spatial(n);
	switch(dir) {
		case YDIR:
			return 1-2*((coord.x)%2);
		case ZDIR:
			return 1-2*((coord.x+coord.y)%2);
		case TDIR:
			return 1-2*((coord.x+coord.y+coord.z)%2);
		default:
			printf("Something bad happened: get_staggered_phase(...) was called without any proper value for dir variable");
	}
}

/** \e Local D_KS working on a particular link (n,t) in a specific direction. The expression of D_KS
 *  for a specific (couple of) site(s) and a specific direction is
 *  \f[
     (D_{KS})_{n,m,\mu}=\frac{1}{2} \eta_\mu(n)\Bigl[U_\mu(n)\,\delta_{n+\hat\mu,m} - U^\dag_\mu(n-\hat\mu)\,\delta_{n-\hat\mu,m}\Bigr]
 *  \f]
 *  and, hence, this function returns the value of the field (D_KS*in) at the site (n,t) [so in the
 *  function only values of the field "in" in the nextneighbour of (n,t) will be needed]: 
 *  \f[
     \bigl[(D_{KS})_\mu\cdot \text{\texttt{in}}\bigr]_n=\frac{1}{2}\eta_\mu(n) \Bigl[U_\mu(n) \cdot\text{\texttt{in}}_{n+\hat\mu} - U^\dag_\mu(n-\hat\mu)\cdot\text{\texttt{in}}_{n-\hat\mu}\Bigr]
 *  \f]
 *
 * @note Again, the staggered phases are included in this function with the help of the function
 *       get_staggered_phase.
 * 
 * @todo Here there are some #ifdef about the chemical potential: they have been copied from the
 *       Wilson code. Therefore they MUST be checked when a chemical potential will be introduced.
 */

//spinor dslash_local_0(__global const spinorfield * const restrict in,__global const ocl_s_gaugefield * const restrict field, const int n, const int t){

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
	U = getSU3(field, get_global_link_pos(dir, n, t));
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
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_acc(out_tmp,chi);
	
	/////////////////////////////////////
	//mu = -0
	//////////////
	nn = get_lower_neighbor_temporal(t);
	plus = get_su3vec_from_field(in, n, nn);
	U = getSU3(field, get_global_link_pos(dir, n, nn));
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
	chi=su3vec_times_complex(chi,bc_tmp);
	out_tmp=su3vec_dim(out_tmp,chi);
	
	///////////////////////////////////
	//multiply by the factor 1/2*eta_t that appears at the beginning of D_KS
	///////////////////////////////////
	eta=0.5*get_staggered_phase(n,dir);
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
	U = getSU3(field, get_global_link_pos(dir, n, t));
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
	U = getSU3(field, get_global_link_pos(dir, nn, t));
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
	U = getSU3(field, get_global_link_pos(dir, n, t));
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
	U = getSU3(field, get_global_link_pos(dir, nn, t));
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
	eta=0.5*get_staggered_phase(n,dir);
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
	U = getSU3(field, get_global_link_pos(dir, n, t));
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
	U = getSU3(field, get_global_link_pos(dir, nn, t));
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
	eta=0.5*get_staggered_phase(n,dir);
	out_tmp = su3vec_times_real(out_tmp, eta); 
	
	return out_tmp;
}

