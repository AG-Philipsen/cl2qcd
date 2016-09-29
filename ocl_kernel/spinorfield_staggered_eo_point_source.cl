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


/*
 * Use the function get_n_eoprec(int spacepos, int timepos) that is in operations_geometry.cl
 * in order to get the eo_superindex to be used in put_su3vec_to_field_eo (see spinorfield_staggered_eo.cl file)
*/



__kernel void create_point_source_stagg_eoprec(__global staggeredStorageType * const restrict inout, int i, int spacepos, int timepos)
{

	const su3vec val;
	val.e0 = 0. ;
	val.e1 = 0. ;
	val.e2 = 0. ;
	
	for(int x = 0; x < NSPACE; x++) 
	{
  		for(int y = 0; y < NSPACE; y++) 
  		{
  			for(int z = 0; z < NSPACE; z++) 
  			{
  				uint3 coord;
			    coord.x = x;
			    coord.y = y;
			    coord.z = z;
			    int iterSpacepos = get_nspace(coord);
			      	
     			for (int iterTimepos = 0; iterTimepos < NTIME_LOCAL; t++)
     			{
			      	uint iterIdx = get_n_eoprec(iterSpacepos, iterTimepos);
      				put_su3vec_to_field_eo(inout, iterIdx, val);
      			}
      		}
      	}
    }

	uint idx = get_n_eoprec(spacepos, timepos);
	hmc_float tmp = 1.;
	switch (i) 
	{
		case 0:
			val.e0 = tmp;
			break;
		case 1:
			val.e1 = tmp;
			break;
		case 2:
			val.e2 = tmp;
			break;
	}
	
	put_su3vec_to_field_eo(inout, idx, val);
	
	return;
}
