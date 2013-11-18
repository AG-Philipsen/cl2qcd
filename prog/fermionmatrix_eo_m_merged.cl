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

__kernel void M_tm_sitediagonal_AND_gamma5_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	//there is a sign change in the factors compared because of the gamma5
	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {-1., 1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);

		putSpinor_eo(out, id_mem, out_tmp);
	}
}

__kernel void M_tm_sitediagonal_minus_AND_gamma5_eo(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	//there is a sign change in the factors compared because of the gamma5
	hmc_complex twistfactor = {1., -mubar_in};
	hmc_complex twistfactor_minus = {-1., -1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);

		putSpinor_eo(out, id_mem, out_tmp);
	}
}
