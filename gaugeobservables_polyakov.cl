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
 * Polyakov calculation kernels
 * NOTE: The reduction used in this kernel is only safe with ls being a power of 2 and bigger than 8!
 */

__kernel void polyakov_reduction(__global hmc_complex* poly_buf,  __global hmc_complex* poly, const uint bufElems)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_complex sum = hmc_complex_zero;
		for (int i = 0; i < bufElems; i++) {
			sum.re += poly_buf[i].re;
			sum.im += poly_buf[i].im;
		}
		*poly = sum;
	}

	return;
}

#if NTIME_GLOBAL == NTIME_LOCAL
/**
 * Calculate the polyakov of the lattice.
 *
 * This method cannot be used in multi-device mode. In that case an alternative approach is required.
 */
__kernel void polyakov(__global Matrixsu3StorageType * field, __global hmc_complex * out, __local hmc_complex * out_loc)
{

	int id;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id_tmp = get_global_id(0);
	int idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_complex tmp_pol;
	hmc_complex tmpcomplex;
	tmp_pol.re = 0.;
	tmp_pol.im = 0.;

	if(idx == 0) {
		out[group_id].re = 0.0f;
		out[group_id].im = 0.0f;
	}

	for(id = id_tmp; id < VOLSPACE; id += global_size) {
		Matrixsu3 prod;
		prod = local_polyakov(field, id);
		tmpcomplex = trace_matrixsu3(prod);
		(tmp_pol).re += tmpcomplex.re / NC;
		(tmp_pol).im += tmpcomplex.im / NC;
	}

	if(local_size == 1) {
		((out))[group_id].re += tmp_pol.re;
		((out))[group_id].im += tmp_pol.im;
	} else {
		//perform reduction
		out_loc[idx].re = tmp_pol.re;
		out_loc[idx].im = tmp_pol.im;
		barrier(CLK_LOCAL_MEM_FENCE);

		//reduction until threads 0-7 hold all partial sums
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 4; cut1 /= 2) {
			for(int i = idx + cut1; i < cut2; i += cut1) {
				((out_loc)[idx]).re +=  ((out_loc)[i]).re;
				((out_loc)[idx]).im +=  ((out_loc)[i]).im;
			}
			barrier(CLK_LOCAL_MEM_FENCE);
			cut2 = cut1;
		}
		//thread 0 sums up the last 8 results and stores them in the global buffer
		if(idx == 0) {
			out[group_id].re = out_loc[0].re + out_loc[1].re + out_loc[2].re + out_loc[3].re +
			  out_loc[4].re + out_loc[5].re + out_loc[6].re + out_loc[7].re;
			out[group_id].im = out_loc[0].im + out_loc[1].im + out_loc[2].im + out_loc[3].im + 
			  out_loc[4].im + out_loc[5].im + out_loc[6].im + out_loc[7].im;
		}
	}
	return;
}

#else

/**
 * Perform the local part of the polyakov calculation.
 */
__kernel void polyakov_md_local(__global Matrixsu3 * local_res_buf, __global const Matrixsu3StorageType * const restrict field)
{
	PARALLEL_FOR(id, VOLSPACE) {
		local_res_buf[id] = local_polyakov(field, id);
	}
}

/**
 * Merge the local polykov results of multiple devices on a single one.
 */
__kernel void polyakov_md_merge(__global hmc_complex * const restrict out, __global const Matrixsu3 * const restrict local_res_bufs, const unsigned num_slices, __local hmc_complex * out_loc)
{
	const size_t local_size = get_local_size(0);
	const size_t idx = get_local_id(0);
	const size_t group_id = get_group_id(0);

	hmc_complex tmp_pol;
	hmc_complex tmpcomplex;
	tmp_pol.re = 0.;
	tmp_pol.im = 0.;

	PARALLEL_FOR(id, VOLSPACE) {
		Matrixsu3 prod = unit_matrixsu3();
		for(unsigned i = 0; i < num_slices; i++) {
			Matrixsu3 tmp = local_res_bufs[i * VOLSPACE + id];
			prod = multiply_matrixsu3(prod, tmp);
		}
		tmpcomplex = trace_matrixsu3(prod);
		tmp_pol.re += tmpcomplex.re / NC;
		tmp_pol.im += tmpcomplex.im / NC;
	}

	//reduction
	if(local_size == 1) {
		out[group_id].re += tmp_pol.re;
		out[group_id].im += tmp_pol.im;
	} else {
		//wait for all threads to end calculations
		barrier(CLK_LOCAL_MEM_FENCE);

		//!!CP: this should be checked by someone else than me
		//perform reduction
		out_loc[idx].re = tmp_pol.re;
		out_loc[idx].im = tmp_pol.im;
		int cut1;
		int cut2 = local_size;
		for(cut1 = local_size / 2; cut1 > 0; cut1 /= 2) {
			for(int i = idx + cut1; i < cut2; i += cut1) {
				((out_loc)[idx]).re +=  ((out_loc)[i]).re;
				((out_loc)[idx]).im +=  ((out_loc)[i]).im;
			}
			//!!CP: is this dangerous inside a for-loop?
			barrier(CLK_LOCAL_MEM_FENCE);
			cut2 = cut1;
		}
		if(idx == 0) {
			out[group_id].re = out_loc[0].re;
			out[group_id].im = out_loc[0].im;
		}
	}
}

#endif
