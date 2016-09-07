
#ifndef _HARDWARE_CODE_UTIL_
#define _HARDWARE_CODE_UTIL_

#include "../../common_header_files/globaldefs.h"
#include "../../common_header_files/types.h"
#include <iostream>
#include <sstream>
#include <string.h>
#include <utility>
//#include "../../host_functionality/host_geometry.h"

namespace hardware {
	namespace code {

	inline size_t getFlopComplexMult() noexcept
	{
		return 6;
	}
	inline size_t getFlopSu3MatrixTimesSu3Matrix() noexcept
	{
		return (getFlopComplexMult() * NC + (NC - 1) * 2) * NC * NC;
	}
	inline size_t getFlopSu3MatrixTrace() noexcept
	{
		return (NC - 1) * 2;
	}
	inline size_t getFlopSu3MatrixTimesSu3Vec() noexcept
	{
		//  1 entry: NC * complex mults and NC-1 complex adds
		//  NC entries total
		return (getFlopComplexMult() * NC + (NC - 1) * 2) * NC;
	}
	inline size_t getFlopSu3VecTimesSu3Vec() noexcept
	{
		// NC * complex_mult + (NC -1) complex adds
		return NC * getFlopComplexMult() + (NC - 1) * 2;
	}
	inline size_t getFlopSu3VecDirectSu3Vec() noexcept
	{
		// NC*NC complex_mult
		return NC * NC * getFlopComplexMult();
	}
	inline size_t getSu3AlgebraSize() noexcept
	{
		return NC * NC - 1;
	}
	inline size_t getFlopSpinorTimesSpinor() noexcept
	{
		//  NDIM * NC * complex_mult + ( NDIM * NC -1 ) complex adds
		return NDIM * NC * getFlopComplexMult() + (NDIM * NC - 1) * 2;
	}
	inline size_t getFlopSpinorSquareNorm() noexcept
	{
		// NDIM * NC * 0.5 complex_mult + ( NDIM * NC -1 ) real adds
		// The 0.5 factor arises from the fact that the result is real and then there is no
		// imaginary part to be calculated.
		return NDIM * NC * getFlopComplexMult() * 0.5 + (NC * NDIM - 1);
	}
	inline size_t getFlopSu3VecSquareNorm() noexcept
	{
		// NDIM * NC * 0.5 complex_mult + ( NDIM * NC -1 ) real adds
		// The 0.5 factor arises from the fact that the result is real and then there is no
		// imaginary part to be calculated.
		return NC * getFlopComplexMult() * 0.5 + (NC - 1);
	}
	}
}
#endif /* _HARDWARE_CODE_UTIL_ */
