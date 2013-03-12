/**
 * Declaration of a utility datatype for fourdimensional sizes
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#endif /* _META_SIZE_4_HPP */
