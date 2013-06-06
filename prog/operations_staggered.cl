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
 
int get_staggered_phase(const int n, const int t, const int dir)
{
	coord_spatial coord;
	coord=get_coord_spatial(n);
	switch(dir) {
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
			printf("Something bad happened: get_staggered_phase(...) was called without any proper value for dir variable");
	}
}