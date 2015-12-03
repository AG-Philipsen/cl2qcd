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

__kernel void create_zslice_source(__global spinor * const restrict b, __global rngStateStorageType * const restrict rngStates, int const zslice)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	spinor out_tmp;
	hmc_complex tmp;
	st_idx pos;
	uint3 coord;

	hmc_float sigma;

	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
	  for(int y = 0; y<NSPACE; y++) {
	  	for(int z = 0; z<NSPACE; z++) {
	      for (int t = 0; t < NTIME_LOCAL; t++){
	      coord.x = id_tmp;
	      coord.y = y;
	      coord.z = z;
	      pos.space = get_nspace(coord);
	      pos.time = t;
	      
	      if( z != zslice )
	      {
	      	out_tmp = set_spinor_zero();
	      }
	      else
	      {
	      //CP: switch between source content
	      switch(SOURCE_CONTENT){
	      case 1:  //"one"
		out_tmp = set_spinor_cold();
		break;
		
	      case 2: //"z4"
		tmp = Z4_complex_number(&rnd);
		out_tmp.e0.e0.re = tmp.re;
		out_tmp.e0.e0.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e0.e1.re = tmp.re;
		out_tmp.e0.e1.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e0.e2.re = tmp.re;
		out_tmp.e0.e2.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e1.e0.re = tmp.re;
		out_tmp.e1.e0.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e1.e1.re = tmp.re;
		out_tmp.e1.e1.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e1.e2.re = tmp.re;
		out_tmp.e1.e2.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e2.e0.re = tmp.re;
		out_tmp.e2.e0.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e2.e1.re = tmp.re;
		out_tmp.e2.e1.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e2.e2.re = tmp.re;
		out_tmp.e2.e2.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e3.e0.re = tmp.re;
		out_tmp.e3.e0.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e3.e1.re = tmp.re;
		out_tmp.e3.e1.im = tmp.im;
		tmp = Z4_complex_number(&rnd);
		out_tmp.e3.e2.re = tmp.re;
		out_tmp.e3.e2.im = tmp.im;
		break;
		
	      case 3: //"gaussian"
		/** @todo what is the norm here? */
		sigma = 0.5;
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
		break;
		
	      default:
		if(id == 0) printf("Problem occured in source kernel: Selected sourcecontent not implemented! Fill with zero...\n");
		out_tmp = set_spinor_zero();
	      }
	      }
	      put_spinor_to_field(out_tmp, b, pos.space, pos.time);
	    }
	  }
	}	    
	}
	prng_storeState(rngStates, &rnd);
	    
	return;
}
	  
