/*
 * Copyright (c) 2011,2012 Christopher Pinke
 * Copyright (c) 2011-2013 Matthias Bach
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file M_tm_minus this is the "total" twisted-mass fermionmatrix (no evenodd) for the lower flavor
 */
__kernel void M_tm_minus(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, __global spinor * const restrict out, hmc_float kappa_in, hmc_float mubar_in)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp, out_tmp2;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {

		/** @todo this must be done more efficient */
		st_index pos = (id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2);

		out_tmp = set_spinor_zero();
		out_tmp2 = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, pos.space, pos.time);
		//Diagonalpart: this is normal tm-diagonal matrix with negative imaginary part
		out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);
		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local_0(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_1(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_2(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_3(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		put_spinor_to_field(out_tmp, out, pos.space, pos.time);
	}
}
