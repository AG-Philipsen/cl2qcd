/*
 * Copyright 2016, Alessandro Sciarra, Tim Breitenfelder
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

/**
 @file fermion-observables
*/

// Correlator is given by:
// C(t)= - 64 * (-1)^t sum_{vec{x}} sum_{c,d} |[D^(-1)_f (vec{x}|0)]_{c,d}|^2
// We drop the factor of (-1)^t since it can be neglected for the statistical evaluation 


// Squarenormcalculation:
//
//hmc_float spinor_staggered_squarenorm(su3vec in)
//
// consider to use this function instead:
// global_squarenorm_staggered_eoprec( __global const staggeredStorageType * const restrict x,
//									  __global hmc_float * const restrict result, __local hmc_float * const restrict result_local )



// this is the pseudoscalar pion correlator in t-direction from pointsources

__kernel void correlator_staggered_ps(__global hmc_float * const restrict out, __global const staggeredStorageType * const restrict phi)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_float summedSquarenorms = 0.;
		uint3 coord;
		int t = id_tmp;
		
		for(coord.x = 0; coord.x < NSPACE; coord.x++) 
		{
			for(coord.y = 0; coord.y < NSPACE; coord.y++) 
			{
				for(coord.z = 0; coord.z < NSPACE; coord.z++) 
				{
					int nspace = get_nspace(coord);
					su3vec tmp = phi[get_n_eoprec(nspace, t)];
					
					// taking squarenorm:
					summedSquarenorms += su3vec_squarenorm(tmp);
				}
			}
		}
		
		//consider taking into account an normalization factor like in wilson case:
		
		//hmc_float fac = NSPACE * NSPACE * NSPACE;
		//this line needs to be modified!
		//out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * correlator / fac;
		
		out[NTIME_OFFSET + id_tmp]= -64 * summedSquarenorms;
		
		
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}
