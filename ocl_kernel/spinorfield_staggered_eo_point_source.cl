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
	int id = get_global_id(0);
	if(id == 0) {
		//LZ: Note that in the kappa-format the source has norm^2 = 1/(2*kappa)
		//we assume everything on the device to be in kappa-normalization
		//going back is achieved by multiplying all fields by sqrt(2*kappa)
		hmc_float tmp = 1.;
		int color     = spinor_color(i);
		int spin      = spinor_spin(i, color);
		int pos       = n;

		spinor site = set_spinor_zero();

		switch (color) {

			case 0:
				switch (spin) {
					case 0:
						site.e0.e0.re = tmp;
						break;
					case 1:
						site.e1.e0.re = tmp;
						break;
					case 2:
						site.e2.e0.re = tmp;
						break;
					case 3:
						site.e3.e0.re = tmp;
						break;
				}
				break;
			case 1:
				switch (spin) {
					case 0:
						site.e0.e1.re = tmp;
						break;
					case 1:
						site.e1.e1.re = tmp;
						break;
					case 2:
						site.e2.e1.re = tmp;
						break;
					case 3:
						site.e3.e1.re = tmp;
						break;
				}
				break;
			case 2:
				switch (spin) {
					case 0:
						site.e0.e2.re = tmp;
						break;
					case 1:
						site.e1.e2.re = tmp;
						break;
					case 2:
						site.e2.e2.re = tmp;
						break;
					case 3:
						site.e3.e2.re = tmp;
						break;
				}
				break;
		}

		putSpinor_eo(out, pos, tmp);
	}
	return;
}
