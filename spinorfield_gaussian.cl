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

__kernel void generate_gaussian_spinorfield(__global spinor * const restrict out, __global rngStateStorageType * const restrict rngStates)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	hmc_complex tmp;
	//sigma has to be 0.5 here
	hmc_float sigma = 0.5;
	spinor out_tmp;

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		/** @todo this must be done more efficient */
		st_index pos = (id_local < SPINORFIELDSIZE_LOCAL / 2) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local - (SPINORFIELDSIZE_LOCAL / 2));

		//CP: there are 12 complex elements in the spinor
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.e0.re = tmp.re;
		out_tmp.e0.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.e1.re = tmp.re;
		out_tmp.e0.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.e2.re = tmp.re;
		out_tmp.e0.e2.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.e0.re = tmp.re;
		out_tmp.e1.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.e1.re = tmp.re;
		out_tmp.e1.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.e2.re = tmp.re;
		out_tmp.e1.e2.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.e0.re = tmp.re;
		out_tmp.e2.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.e1.re = tmp.re;
		out_tmp.e2.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.e2.re = tmp.re;
		out_tmp.e2.e2.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e3.e0.re = tmp.re;
		out_tmp.e3.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e3.e1.re = tmp.re;
		out_tmp.e3.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e3.e2.re = tmp.re;
		out_tmp.e3.e2.im = tmp.im;

		//multiply by sigma
		out_tmp = real_multiply_spinor(out_tmp, sqrt(sigma));

		put_spinor_to_field(out_tmp, out, pos.space, pos.time);
	}

	prng_storeState(rngStates, &rnd);
}
