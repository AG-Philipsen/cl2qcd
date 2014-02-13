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

__kernel void create_volume_source_stagg_eoprec(__global staggeredStorageType * const restrict inout, __global rngStateStorageType * const restrict rngStates)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);
	su3vec out_tmp;
	hmc_complex tmp;
	
	for(int n = id; n < EOPREC_SPINORFIELDSIZE_MEM; n += global_size) {
	  switch(SOURCE_CONTENT){
		case 1:  //"one"
		  out_tmp = set_su3vec_cold();
		  break;
		case 2: //"z4"
		  tmp = Z4_complex_number(&rnd);
		  out_tmp.e0.re = tmp.re;
		  out_tmp.e0.im = tmp.im;
		  tmp = Z4_complex_number(&rnd);
		  out_tmp.e0.re = tmp.re;
		  out_tmp.e0.im = tmp.im;
		  tmp = Z4_complex_number(&rnd);
		  out_tmp.e0.re = tmp.re;
		  out_tmp.e0.im = tmp.im;
		  break;
		case 3: //"gaussian"
		  /** @todo what is the norm here? */
		  tmp = gaussianNormalPair(&rnd);
		  out_tmp.e0.re = tmp.re;
		  out_tmp.e0.im = tmp.im;
		  tmp = gaussianNormalPair(&rnd);
		  out_tmp.e0.re = tmp.re;
		  out_tmp.e0.im = tmp.im;
		  tmp = gaussianNormalPair(&rnd);
		  out_tmp.e0.re = tmp.re;
		  out_tmp.e0.im = tmp.im;
		  //multiply by sigma = 0.5f
		  out_tmp = su3vec_times_real(out_tmp, sqrt(0.5f));
		  break;
		default:
		  if(id == 0) printf("Problem occured in source kernel: Selected sourcecontent not implemented! Fill with zero...\n");
		  out_tmp = set_su3vec_zero();
	  }
	  put_su3vec_to_field_eo(inout, n, out_tmp);
	}
	prng_storeState(rngStates, &rnd);
}



