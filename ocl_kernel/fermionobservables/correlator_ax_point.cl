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

/**
 @file fermion-observables
*/

hmc_complex calculate_ax_correlator(const spinor in1, const spinor in2, const spinor in3, const spinor in4)
{
	hmc_complex restmp;
	hmc_complex correlator;
	correlator.re = 0.0f;
	correlator.im = 0.0f;

	restmp = su3vec_scalarproduct(in1.e0, in2.e1);
	correlator.re += restmp.re;
	correlator.im += restmp.im;

	restmp = su3vec_scalarproduct(in1.e1, in2.e0);
	correlator.re += restmp.re;
	correlator.im += restmp.im;

	restmp = su3vec_scalarproduct(in1.e2, in2.e3);
	correlator.re -= restmp.re;
	correlator.im -= restmp.im;

	restmp = su3vec_scalarproduct(in1.e3, in2.e2);
	correlator.re -= restmp.re;
	correlator.im -= restmp.im;

	restmp = su3vec_scalarproduct(in3.e0, in4.e1);
	correlator.re -= restmp.re;
	correlator.im -= restmp.im;

	restmp = su3vec_scalarproduct(in3.e1, in4.e0);
	correlator.re -= restmp.re;
	correlator.im -= restmp.im;

	restmp = su3vec_scalarproduct(in3.e2, in4.e3);
	correlator.re += restmp.re;
	correlator.im += restmp.im;

	restmp = su3vec_scalarproduct(in3.e3, in4.e2);
	correlator.re += restmp.re;
	correlator.im += restmp.im;
	
	return correlator;
}

// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma1)
__kernel void correlator_ax_z(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict phi2, __global const spinor * const restrict phi3, __global const spinor * const restrict phi4)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME_LOCAL; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					hmc_complex cortmp;
					spinor tmp1 = phi1[get_pos(nspace, t)];
					spinor tmp2 = phi2[get_pos(nspace, t)];
					spinor tmp3 = phi3[get_pos(nspace, t)];
					spinor tmp4 = phi4[get_pos(nspace, t)];
					
					cortmp = calculate_ax_correlator(tmp1, tmp2, tmp3, tmp4);
					correlator.re += cortmp.re;
					correlator.im += cortmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME_GLOBAL;
		out[id_tmp] += - 2. * KAPPA * 2.*KAPPA * 2.*correlator.re / fac;	//why the minus sign?
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma1)
__kernel void correlator_ax_t(__global hmc_float * const restrict out, __global const spinor * const restrict phi1, __global const spinor * const restrict phi2, __global const spinor * const restrict phi3, __global const spinor * const restrict phi4)
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
					spinor tmp1 = phi1[get_pos(nspace, t)];
					spinor tmp2 = phi2[get_pos(nspace, t)];
					spinor tmp3 = phi3[get_pos(nspace, t)];
					spinor tmp4 = phi4[get_pos(nspace, t)];
					
					cortmp = calculate_ax_correlator(tmp1, tmp2, tmp3, tmp4);
					correlator.re += cortmp.re;
					correlator.im += cortmp.im;
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[NTIME_OFFSET + id_tmp] += - 2. * KAPPA * 2.*KAPPA * 2. * correlator.re / fac;	//why the minus sign?
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}