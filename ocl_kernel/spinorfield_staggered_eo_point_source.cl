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


__kernel void create_point_source_stagg_eoprec(__global staggeredStorageType * const restrict inout, int i, int spacepos, int timepos)
{
	int id = get_global_id(0);
    int global_size = get_global_size(0);
    
    for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
        put_su3vec_to_field_eo(inout, id_mem, set_su3vec_zero());
        uint idx = get_n_eoprec(spacepos, timepos);
        //The following if is to avoid race conditions. The source is on a site and the thread that set to zero
        //that site must be the same that set it to non-zero later. If two different thread do this, it could
        //be that the thread that set it to non-zero acts before the thread that set it to zero!
    	if(id_mem == idx) {
    	    //No default case in switch since it is checked before enqueuing the kernel that i=0,1,2 
    	    su3vec val = set_su3vec_zero();
    	    hmc_float tmp = 1.;
    	    switch (i) 
    	    {
    	        case 0:
    	            val.e0.re = tmp;
    	            break;
    	        case 1:
    	            val.e1.re = tmp;
    	            break;
    	        case 2:
    	            val.e2.re = tmp;
    	            break;
    	    }
    	    put_su3vec_to_field_eo(inout, idx, val);
    	}
    }

	return;
}
