/** @file
 * General utilities for staggered code (e.g. calculation of staggered phases)
 */

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
         \f]
 *       where is CRUCIAL to emphasize that
 *       \f[
             \begin{aligned}
               \mu=1 &\quad\mbox{means}\quad x=n_1 \\
               \mu=2 &\quad\mbox{means}\quad y=n_2 \\
               \mu=3 &\quad\mbox{means}\quad z=n_3 \\
               \mu=4 &\quad\mbox{means}\quad t=n_4 
             \end{aligned}
         \f]
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
 * @param n this is the spatial superindex of the site (see st_idx.space in operation_geometry.cl)
 * @param dir this is the direction of the staggered phase. To be automatically coherent with
 *             the choice of labels made in operation_geometry.cl, it can be YDIR, ZDIR or TDIR.
 *      
 */
int get_staggered_phase(const int n, const int dir)
{
	coord_spatial coord = get_coord_spatial(n);
	switch(dir) {
// 		case XDIR:
// 			return 1;
		case YDIR:
			return 1-2*((coord.x)%2);
		case ZDIR:
			return 1-2*((coord.x+coord.y)%2);
		case TDIR:
			//if(t!=(NTIME-1))
			return 1-2*((coord.x+coord.y+coord.z)%2);
			//else
			//return -(1-2*((coord.x+coord.y+coord.z)%2));
		default:
			return 1;
	}
}

/**
 * This function returns the staggered phase modified in order to include in it
 * the boundary conditions, as in the Wilson code. This means that all staggered
 * phases are multiplied by exp(i*theta_mu/N_mu), where N_mu is the lattice extension
 * in the mu direction. For this reason the return value is a complex number.
 */
hmc_complex get_modified_stagg_phase(const int n, const int dir)
{
	coord_spatial coord = get_coord_spatial(n);
	hmc_complex out;
	switch(dir) { //Note that the 2. (instead of 2) is crucial to have a float
	  case YDIR:
	    out.re = (1-2.*((coord.x)%2))*SPATIAL_RE;
	    out.im = (1-2.*((coord.x)%2))*SPATIAL_IM;
	    break;
	  case ZDIR:
	    out.re = (1-2.*((coord.x+coord.y)%2))*SPATIAL_RE;
	    out.im = (1-2.*((coord.x+coord.y)%2))*SPATIAL_IM;
	    break;
	  case TDIR:
	    out.re = (1-2.*((coord.x+coord.y+coord.z)%2))*TEMPORAL_RE;
	    out.im = (1-2.*((coord.x+coord.y+coord.z)%2))*TEMPORAL_IM;
	    break;
	  default:
	    printf("Assuming dir=XDIR");
	    out.re = 1.*SPATIAL_RE;
	    out.im = 1.*SPATIAL_IM;
	    break;
	}

	return out;
}


/**
 * This function returns the staggered phase modified in order to include in it
 * the boundary conditions. This means that the staggered phases at the last site
 * in the mu-dir must be multiplied by exp(i*theta_mu). For this reason the return
 * value is a complex number and we have also to pass to the function the time coordinate
 * in order to judge whether we are on the last site in the temporal direction.
 */
/*
//In case the following function is uncommented, these options have to be passed
//to the kernels that use the function.
//These 4 parameters are needed to modify staggered phases and then to impose BC
//  options << " -D COS_THETAS=" << cos(params.get_theta_fermion_spatial() * PI);
//  options << " -D SIN_THETAS=" << sin(params.get_theta_fermion_spatial() * PI);
//  options << " -D COS_THETAT=" << cos(params.get_theta_fermion_temporal() * PI);
//  options << " -D SIN_THETAT=" << sin(params.get_theta_fermion_temporal() * PI);

hmc_complex get_mod_stagg_phase(const int n, const int t, const int dir)
{
	int ph;
	coord_spatial coord = get_coord_spatial(n);
	hmc_complex out;
	switch(dir) {
	  case YDIR:
		ph = 1-2*((coord.x)%2);
		if(coord.y == (NSPACE-1)) {
			out.re = ph * COS_THETAS;
			out.im = ph * SIN_THETAS;
		} else {
			out.re = ph;
			out.im = 0.;
		}
		break;
	  case ZDIR:
		ph = 1-2*((coord.x+coord.y)%2);
		if(coord.z == (NSPACE-1)) {
			out.re = ph * COS_THETAS;
			out.im = ph * SIN_THETAS;
		} else {
			out.re = ph;
			out.im = 0.;
		}
		break;
	  case TDIR:
		ph = 1-2*((coord.x+coord.y+coord.z)%2);
		if(t == (NTIME_GLOBAL-1)) {
			out.re = ph * COS_THETAT;
			out.im = ph * SIN_THETAT;
		} else {
			out.re = ph;
			out.im = 0.;
		}
		break;
	  default:
		printf("Assuming dir=XDIR");
		if(coord.x == (NSPACE-1)) {
			out.re = COS_THETAS;
			out.im = SIN_THETAS;
		} else {
			out.re = 1.;
			out.im = 0.;
		}
		break;
	}
	return out;
}
*/




