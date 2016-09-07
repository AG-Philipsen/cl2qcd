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

__kernel void create_point_source(__global spinor * const restrict b, int i, int spacepos, int timepos)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);
	
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
	  for(int y = 0; y<NSPACE; y++) {
	  	for(int z = 0; z<NSPACE; z++) {
	      for (int t = 0; t < NTIME_LOCAL; t++){
			uint3 coord;
	      	coord.x = id_tmp;
	      	coord.y = y;
	      	coord.z = z;
	      	int posSpace = get_nspace(coord);
	      
			put_spinor_to_field(set_spinor_zero(), b, posSpace, t);  
	      }
	    }
	  }
	}

	if(id == 0) {
		if (SOURCE_CONTENT != 1) //"1" is "one" 
		{
			printf("Problem occured in source kernel: Selected sourcecontent not implemented! Fill with zero...\n");
			return;
		}
	
		//LZ: note that the conversion from m to kappa works as
		//   M(m) * phi = eta <=> 1/(2kappa) * M(k) * phi = eta
		//thus we can keep everything as in the orginal basis and only have
		//to multiply the resulting field phi by 2kappa
		hmc_float tmp = 1.;
		int color = spinor_color(i);
		int spin = spinor_spin(i, color);
		int pos = get_pos(spacepos, timepos);
		b[pos] = set_spinor_zero();
		switch (color) {

			case 0:
				switch (spin) {
					case 0:
						(b[pos].e0).e0.re = tmp;
						break;
					case 1:
						(b[pos].e1).e0.re = tmp;
						break;
					case 2:
						(b[pos].e2).e0.re = tmp;
						break;
					case 3:
						(b[pos].e3).e0.re = tmp;
						break;
				}
				break;
			case 1:
				switch (spin) {
					case 0:
						(b[pos].e0).e1.re = tmp;
						break;
					case 1:
						(b[pos].e1).e1.re = tmp;
						break;
					case 2:
						(b[pos].e2).e1.re = tmp;
						break;
					case 3:
						(b[pos].e3).e1.re = tmp;
						break;
				}
				break;
			case 2:
				switch (spin) {
					case 0:
						(b[pos].e0).e2.re = tmp;
						break;
					case 1:
						(b[pos].e1).e2.re = tmp;
						break;
					case 2:
						(b[pos].e2).e2.re = tmp;
						break;
					case 3:
						(b[pos].e3).e2.re = tmp;
						break;
				}
				break;
		}
	}
	return;
}
