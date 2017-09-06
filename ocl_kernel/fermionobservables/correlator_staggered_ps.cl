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
// We drop the prefactor of -64*(-1)^t in the following since it can be neglected for the final mass extraction 

__kernel void correlator_staggered_ps(__global hmc_float * const restrict correlator, __global const staggeredStorageType * const restrict invertedSourceEven,
 																			   		  __global const staggeredStorageType * const restrict invertedSourceOdd) 
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	for(int id_tmp = id; id_tmp < NTIME_LOCAL; id_tmp += global_size ) 
	{
		hmc_float summedSquarenorms = 0.;
		uint3 coord;
		
		for(coord.x = 0; coord.x < NSPACE; coord.x++ ) 
		{
			for(coord.y = 0; coord.y < NSPACE; coord.y++ ) 
			{
				for(coord.z = 0; coord.z < NSPACE; coord.z++ ) 
				{
					int nspace = get_nspace(coord);
					int tmp_idx = get_n_eoprec(nspace, id_tmp);				
					const bool sourceOnEvenSite = ((coord.x+coord.y+coord.z+id_tmp)%2 == 0) ? true : false;
					su3vec temporalField;					
					
					if(sourceOnEvenSite)
					{
						temporalField = get_su3vec_from_field_eo(invertedSourceEven, tmp_idx);
					}
					else
					{
						temporalField = get_su3vec_from_field_eo(invertedSourceOdd, tmp_idx);
					}
					
					summedSquarenorms += su3vec_squarenorm(temporalField);
				}
			}
		}

		hmc_float spatialVolume = NSPACE * NSPACE * NSPACE;
		correlator[NTIME_OFFSET + id_tmp] += summedSquarenorms / spatialVolume ;		
	}
}