/** @file
 * Operators und utility functions for the custom types
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "type_ops.hpp"
#include <cstdlib>

template<> void fill(hmc_complex* array, size_t num_elems, int seed)
{
	srand48(seed);
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {drand48(), drand48()};
	}
}

template<> void fill(Matrixsu3* array, size_t num_elems, int seed)
{
	srand48(seed);
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() },
			{ drand48(), drand48() }
		};
	}
}

template<> void fill(spinor* array, size_t num_elems, int seed)
{
	srand48(seed);
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
			{
				{ drand48(), drand48() },
				{ drand48(), drand48() },
				{ drand48(), drand48() },
			}, {
				{ drand48(), drand48() },
				{ drand48(), drand48() },
				{ drand48(), drand48() },
			}, {
				{ drand48(), drand48() },
				{ drand48(), drand48() },
				{ drand48(), drand48() }
			}, {
				{ drand48(), drand48() },
				{ drand48(), drand48() },
				{ drand48(), drand48() }
			}
		};
	}
}

template<> void fill(su3vec* array, size_t num_elems, int seed)
{
	srand48(seed);
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
				{ drand48(), drand48() },
				{ drand48(), drand48() },
				{ drand48(), drand48() },
		};
	}
}

template<> void fill(ae* array, size_t num_elems, int seed)
{
	srand48(seed);
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
			drand48(),
			drand48(),
			drand48(),
			drand48(),
			drand48(),
			drand48(),
			drand48(),
			drand48()
		};
	}
}

template<> size_t get_flops<hmc_complex, complexconj>()
{
	return 0;
}
template<> size_t get_flops<hmc_complex, complexmult>()
{
	return 6;
}
template<> size_t get_flops<hmc_complex, complexadd>()
{
	return 2;
}
template<> size_t get_flops<hmc_complex, complexsubtract>()
{
	return 2;
}
template<> size_t get_flops<hmc_complex, complexdivide>()
{
	return 11;
}

template<> size_t get_read_write_size<hmc_complex, complexconj>()
{
	return 4*8;
}
template<> size_t get_read_write_size<hmc_complex, complexmult>()
{
	return 6*8;
}
template<> size_t get_read_write_size<hmc_complex, complexadd>()
{
	return 6*8;
}
template<> size_t get_read_write_size<hmc_complex, complexsubtract>()
{
	return 6*8;
}
template<> size_t get_read_write_size<hmc_complex, complexdivide>()
{
	return 6*8;
}
