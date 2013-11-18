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
 * Device code for operations on the spinor field
 */

spinor get_spinor_from_field(__global const spinor * const restrict in, const int n, const int t)
{
	int pos = get_pos(n, t);
	spinor out;
	out = in[pos];
	return out;
}

void put_spinor_to_field(const spinor in, __global spinor * const restrict out, const int n, const int t)
{
	int pos = get_pos(n, t);
	out[pos] = in;
}
