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

//Here, two cases are possible:
//evenodd = ODD or EVEN
//	ODD corresponds to the D_oe case: Dslash acts on even indices (the "x+mu" in the formulae) and the function
//	saves the outcoming spinor with an odd index.
//	EVEN is then D_eo.

void dslash_eo_for_site(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in, st_idx const pos)
{
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

	putSpinor_eo(out, get_eo_site_idx_from_st_idx(pos), out_tmp);
}

__kernel void dslash_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	PARALLEL_FOR(id_local, EOPREC_SPINORFIELDSIZE_LOCAL) {
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		dslash_eo_for_site(in, out, field, evenodd, kappa_in, pos);
	}
}

#define REQD_HALO_WIDTH 1
#define HALO_VOL (VOLSPACE / 2 * REQD_HALO_WIDTH)

__kernel void dslash_eo_inner(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	size_t id_local;
	PARALLEL_FOR(id_loop, EOPREC_SPINORFIELDSIZE_LOCAL - (2 * HALO_VOL)) {
		// note that the scheme we are generating positions will no longer work for spatial seperation!
		id_local = id_loop + HALO_VOL; // boost position by halo width
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		dslash_eo_for_site(in, out, field, evenodd, kappa_in, pos);
	}
}

__kernel void dslash_eo_boundary(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, __global const Matrixsu3StorageType * const restrict field, const int evenodd, hmc_float kappa_in)
{
	size_t id_local;
	PARALLEL_FOR(id_loop, 2 * HALO_VOL) {
		id_local = (id_loop < HALO_VOL) ? id_loop : (EOPREC_SPINORFIELDSIZE_LOCAL - HALO_VOL + (id_loop - HALO_VOL));
		// note that the scheme we are generating positions will no longer work for spatial seperation!
		st_idx pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		dslash_eo_for_site(in, out, field, evenodd, kappa_in, pos);
	}
}
