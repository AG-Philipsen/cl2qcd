/** @file
 * Operators for the custom types
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _META_TYPE_OPS_
#define _META_TYPE_OPS_

#include "../types.h"

#include <ostream>

inline bool operator==(const hmc_complex& left, const hmc_complex& right)
{
	return (left.re == right.re && left.im == right.im);
}
inline bool operator!=(const hmc_complex& left, const hmc_complex& right)
{
	return (left.re != right.re || left.im != right.im);
}
std::ostream& operator<<(std::ostream& os, const hmc_complex& data)
{
	return os << '(' << data.re << ',' << data.im << "i)";
}

inline bool operator==(const Matrixsu3& left, const Matrixsu3& right)
{
	return (
			left.e00 == right.e00 && left.e01 == right.e01 && left.e02 == right.e02
		 && left.e10 == right.e10 && left.e11 == right.e11 && left.e12 == right.e12
		 && left.e20 == right.e20 && left.e21 == right.e21 && left.e22 == right.e22
	);
}
inline bool operator!=(const Matrixsu3& left, const Matrixsu3& right)
{
	return (
			left.e00 != right.e00 || left.e01 != right.e01 || left.e02 != right.e02
		 || left.e10 != right.e10 || left.e11 != right.e11 || left.e12 != right.e12
		 || left.e20 != right.e20 || left.e21 != right.e21 || left.e22 != right.e22
	);
}
std::ostream& operator<<(std::ostream& os, const Matrixsu3& data)
{
	return os << '{' << data.e00 << ',' << data.e01 << ',' << data.e02 << ';'
	                 << data.e10 << ',' << data.e11 << ',' << data.e12 << ';'
	                 << data.e20 << ',' << data.e21 << ',' << data.e22 << '}';
}


#endif /* _META_TYPE_OPS */
