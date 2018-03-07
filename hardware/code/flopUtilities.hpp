/*
 * Copyright (c) 2015 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HARDWARE_CODE_UTIL_
#define _HARDWARE_CODE_UTIL_

#include "../../common_header_files/globaldefs.hpp"
#include "../../common_header_files/types.hpp"
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
