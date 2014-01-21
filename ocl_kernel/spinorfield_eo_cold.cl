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

__kernel void set_eoprec_spinorfield_cold(__global spinorStorageType * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	//NOTE: the normalization is the same as in the above case, putting |phi|^2 of the WHOLE field to 1!!
	hmc_float norm = 1. / sqrt((hmc_float)(12. * VOL4D_GLOBAL));
	for(int n = id; n < EOPREC_SPINORFIELDSIZE_MEM; n += global_size) {
		putSpinor_eo(out, n, real_multiply_spinor(set_spinor_cold(), norm));
	}
}
