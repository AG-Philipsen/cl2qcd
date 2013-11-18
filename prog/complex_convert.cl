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
 * Conversion of complex numbers
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

/**
 * Converts a float to a complex number
 */
__kernel void convert_float_to_complex(__global hmc_float const * const restrict in, __global hmc_complex * const restrict out)
{
	if(get_global_id(0) == 0) {
		hmc_complex tmp;
		tmp.re = *in;
		tmp.im = 0.f;
		*out = tmp;
	}
}

