/**
 * Declaration of a utility datatype for fourdimensional sizes
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _META_SIZE_4_HPP
#define _META_SIZE_4_HPP

#include <ostream>

struct size_4 {
	unsigned x, y, z, t;

	size_4(unsigned x, unsigned y, unsigned z, unsigned t) : x(x), y(y), z(z), t(t) { };

};

inline std::ostream& operator<<(std::ostream& os, const size_4& data)
{
	return os << '(' << data.x << ", " << data.y << ", " << data.z << ", " << data.t << ')';
}

inline unsigned get_vol4d(size_4 data) {
	return data.x * data.y * data.z * data.t;
}

#endif /* _META_SIZE_4_HPP */
