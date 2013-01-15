/** @file
 * Operators und utility functions for the custom types
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _META_TYPE_OPS_
#define _META_TYPE_OPS_

#include "../types.h"
#include "../types_fermions.h"
#include "../operations_complex.h"

#include <ostream>

inline bool operator==(const hmc_complex& left, const hmc_complex& right)
{
	return (left.re == right.re && left.im == right.im);
}
inline bool operator!=(const hmc_complex& left, const hmc_complex& right)
{
	return (left.re != right.re || left.im != right.im);
}
inline std::ostream& operator<<(std::ostream& os, const hmc_complex& data)
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
inline std::ostream& operator<<(std::ostream& os, const Matrixsu3& data)
{
	return os << '{' << data.e00 << ',' << data.e01 << ',' << data.e02 << ';'
	                 << data.e10 << ',' << data.e11 << ',' << data.e12 << ';'
	                 << data.e20 << ',' << data.e21 << ',' << data.e22 << '}';
}

inline bool operator==(const su3vec& left, const su3vec& right)
{
	return (
			left.e0 == right.e0 && left.e1 == right.e1 && left.e2 == right.e2
	);
}
inline bool operator!=(const su3vec& left, const su3vec& right)
{
	return (
			left.e0 != right.e0 || left.e1 != right.e1 || left.e2 != right.e2
	);
}
inline std::ostream& operator<<(std::ostream& os, const su3vec& data)
{
	return os << '{' << data.e0 << ',' << data.e1 << ',' << data.e2 << '}';
}

inline bool operator==(const spinor& left, const spinor& right)
{
	return (
			left.e0 == right.e0 && left.e1 == right.e1 && left.e2 == right.e2 && left.e3 == right.e3
	);
}
inline bool operator!=(const spinor& left, const spinor& right)
{
	return (
			left.e0 != right.e0 || left.e1 != right.e1 || left.e2 != right.e2 || left.e3 != right.e3
	);
}
inline std::ostream& operator<<(std::ostream& os, const spinor& data)
{
	return os << '{' << data.e0 << ',' << data.e1 << ',' << data.e2 << ',' << data.e3 << '}';
}

inline bool operator==(const ae& left, const ae& right)
{
	return (
			left.e0 == right.e0 && left.e1 == right.e1 && left.e2 == right.e2 && left.e3 == right.e3 && left.e4 == right.e4 && left.e5 == right.e5 && left.e6 == right.e6 && left.e7 == right.e7
	);
}
inline bool operator!=(const ae& left, const ae& right)
{
	return (
			left.e0 != right.e0 || left.e1 != right.e1 || left.e2 != right.e2 || left.e3 != right.e3 || left.e4 != right.e4 || left.e5 != right.e5 || left.e6 != right.e6 || left.e7 != right.e7
	);
}
inline std::ostream& operator<<(std::ostream& os, const ae& data)
{
	return os << '{' << data.e0 << ',' << data.e1 << ',' << data.e2 << ',' << data.e3 << ',' << data.e4 << ',' << data.e5 << ',' << data.e7 << '}';
}

template<typename T> inline void fill(T* array, size_t num_elems, int seed = 0)
{
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = i + seed;
	}
}
template<> void fill(hmc_complex* array, size_t num_elems, int seed);
template<> void fill(Matrixsu3* array, size_t num_elems, int seed);
template<> void fill(spinor* array, size_t num_elems, int seed);
template<> void fill(ae* array, size_t num_elems, int seed);

/*
 * OP counts for complex operations
 */
template<typename S, S (*T)(S)> size_t get_flops();
template<typename S, S (*T)(S,S)> size_t get_flops();
template<> size_t get_flops<hmc_complex, complexconj>();
template<> size_t get_flops<hmc_complex, complexmult>();
template<> size_t get_flops<hmc_complex, complexadd>();
template<> size_t get_flops<hmc_complex, complexsubtract>();
template<> size_t get_flops<hmc_complex, complexdivide>();

#endif /* _META_TYPE_OPS */
