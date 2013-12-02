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
 * Device code for operations on the spinor_staggered field
 */

/**
 * This functions returns the value of the field "in" at the site (n,t)
 */
su3vec get_su3vec_from_field(__global const su3vec * const restrict in, const int n, const int t)
{
	int pos = get_pos(n, t);
	su3vec out;
	out = in[pos];
	return out;
}

/**
 * This functions uploads the value of the field "out" at the site (n,t) with the value "in"
 */
void put_su3vec_to_field(const su3vec in, __global su3vec * const restrict out, const int n, const int t)
{
	int pos = get_pos(n, t);
	out[pos] = in;
}
