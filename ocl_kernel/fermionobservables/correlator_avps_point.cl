/*
 * Copyright 2015 Paul Frederik Depta
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

hmc_complex calculate_avps_correlator(const spinor in)
{
	hmc_complex restmp;
	hmc_complex correlator;
	correlator.re = 0.0f;
	correlator.im = 0.0f;
	
	restmp = su3vec_scalarproduct(in.e0, in.e2);
	correlator.re -= restmp.re;
	correlator.im -= restmp.im;

	restmp = su3vec_scalarproduct(in.e1, in.e3);
	correlator.re -= restmp.re;
	correlator.im -= restmp.im;
		
	return correlator;
}

/*	this is the axial vector 4th component with pseudoscalar (AVPS) correlator in t-direction 

	C = tr((D^(-1)(n|m))^dagger * gamma4 * D^(-1)(n|m)), where m is the source position and n is the sink position
	
                  |0   0  -1   0|
	gamma4 =  |0   0   0  -1|
                  |-1  0   0   0|
                  |0  -1   0   0|

	D^(-1)(n|m) = (phi) with summation over second spin-color index elsewhere.
*/ 
__kernel void correlator_avps_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					hmc_complex cortmp;
					spinor tmp = phi[get_pos(nspace, t)];
					
					cortmp = calculate_avps_correlator(tmp);
					correlator.re += cortmp.re;
					correlator.im += cortmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += 2. * KAPPA * 2. * KAPPA * 2. * correlator.re / fac;
	}
}