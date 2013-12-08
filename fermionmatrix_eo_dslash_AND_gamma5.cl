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

//this is the kernel merging dslash_eo and gamma5

//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.

__kernel void dslash_AND_gamma5_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);

		spinor out_tmp = set_spinor_zero();
		spinor out_tmp2;

		//calc dslash (this includes mutliplication with kappa)

		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, TDIR, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, XDIR, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, YDIR, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_eoprec_unified_local(in, field, pos, ZDIR, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		//gamma5 part
		out_tmp = gamma5_local(out_tmp);

		putSpinor_eo(out, get_eo_site_idx_from_st_idx(pos), out_tmp);
	}
}
