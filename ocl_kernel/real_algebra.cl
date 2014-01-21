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
 * Device code implementing real numbers algebra functionalities
 */

__kernel void update_zeta_cgm(__global const hmc_float * zeta_prev, __global const hmc_float * zeta_prev_prev, __global const hmc_float * sbeta_prev, __global const hmc_float * sbeta_pres, __global const hmc_float * salpha_prev, __global const hmc_float * sigma, const int numeq, __global hmc_float * out)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	
	for(int id_mem = id; id_mem < numeq; id_mem += global_size) 
	   out[id_mem] = update_zeta_cgm_alg(zeta_prev[id_mem], zeta_prev_prev[id_mem], *sbeta_prev,
					      *sbeta_pres, *salpha_prev, sigma[id_mem]);
	
}


__kernel void update_beta_cgm(__global const hmc_float * sbeta_pres, __global const hmc_float * zeta_pres, __global const hmc_float * zeta_prev, const int numeq, __global hmc_float * out)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	
	for(int id_mem = id; id_mem < numeq; id_mem += global_size) 
	   out[id_mem] = update_beta_cgm_alg(*sbeta_pres, zeta_pres[id_mem], zeta_prev[id_mem]);
	
}


__kernel void update_alpha_cgm(__global const hmc_float * salpha_pres, __global const hmc_float * zeta_pres, __global const hmc_float * beta_pres, __global const hmc_float * zeta_prev, __global const hmc_float * sbeta_pres, const int numeq, __global hmc_float * out)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	
	for(int id_mem = id; id_mem < numeq; id_mem += global_size) 
	   out[id_mem] = update_alpha_cgm_alg(*salpha_pres, zeta_pres[id_mem], beta_pres[id_mem],
					        zeta_prev[id_mem], *sbeta_pres);
	
}




