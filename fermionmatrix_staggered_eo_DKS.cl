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
 * Even-odd or Odd-even block of the staggered standard Dirac operator
 *
 * Following the standard notation (see file fermionmatrix_staggered_M.cl for details)
 * we can easily see that M = m + D_KS is exactly the sum of a part that connects sites
 * of the same parity (i.e. m that does not contain links), and a part that connects sites
 * of opposite parity (i.e. D_KS). The latter can be splitted in the two operators Doe and
 * Deo that connects respectively even sites to odd sites and odd sites to even sites.
 * For instance Doe acting on an even field returns an odd one and vice-versa.
 * 
 * Here in this kernel we make the Dxy operator act on the field "in" getting the result in the
 * field "out". I wrote on purpose Dxy meaning that this kernel can be both Doe and Deo. To switch
 * between the two one must use the variable evenodd that ALWAYS corresponds to the parity of the
 * OUTPUT field. 
 * 
 * @param in The even/odd input staggered field
 * @param out The odd/even output staggered field
 * @param field The gauge configuration
 * @param evenodd Parameter to switch between Doe (evenodd==ODD) and Deo (evenodd==EVEN)
 * 
 * @note This kernel takes into account the staggered phases eta_mu(n), since they
 *       will be NOT included in links (as in the phi-algorithm from Gottlieb and Toussaint).
 *       The reason why we make such a choice is that we can preserve gauge configurations
 *       (and all related code) as it is, without changing the links in the staggered algorithm.
 * \par
 * @note In this kernel \em even-odd \em preconditioning is assumed: vectors "in" and "out" have
 *       VOL4D/2 components.
 */

__kernel void D_KS_eo(__global const staggeredStorageType * const restrict in, __global staggeredStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd)
{
	PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
		st_idx pos = (evenodd == EVEN) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);

		su3vec out_tmp = set_su3vec_zero();
		su3vec out_tmp2;

		//Non-diagonal part: calc D_KS (here if it is Doe or Deo is automatic
		//                              thanks to the if above to set pos)
		for(dir_idx dir = 0; dir < 4; ++dir) {
			out_tmp2 = D_KS_eo_local(in, field, pos, dir);
			out_tmp = su3vec_acc(out_tmp, out_tmp2);
		}


		put_su3vec_to_field_eo(out, get_eo_site_idx_from_st_idx(pos), out_tmp);
	}
}
