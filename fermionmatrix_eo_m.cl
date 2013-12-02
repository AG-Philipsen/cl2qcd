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

__kernel void M_tm_sitediagonal(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
		putSpinor_eo(out, id_mem, out_tmp);
	}
}

__kernel void M_tm_inverse_sitediagonal(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	spinor out_tmp;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		//Diagonalpart, here the twisted factor give the inverse matrix:
		out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);
		hmc_float denom = 1. / (1. + mubar_in * mubar_in);
		out_tmp = real_multiply_spinor(out_tmp, denom);
		putSpinor_eo(out, id_mem, out_tmp);
	}
}

__kernel void M_tm_sitediagonal_minus(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);
		putSpinor_eo(out, id_mem, out_tmp);
	}
}

__kernel void M_tm_inverse_sitediagonal_minus(__global const spinorStorageType * const restrict in, __global spinorStorageType * const restrict out, hmc_float mubar_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	spinor out_tmp;
	spinor plus;
	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = getSpinor_eo(in, id_mem);
		//Diagonalpart, here the twisted factor give the inverse matrix:
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
		hmc_float denom = 1. / (1. + mubar_in * mubar_in);
		out_tmp = real_multiply_spinor(out_tmp, denom);
		putSpinor_eo(out, id_mem, out_tmp);
	}
}

