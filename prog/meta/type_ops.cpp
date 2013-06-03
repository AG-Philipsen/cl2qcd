/** @file
 * Operators und utility functions for the custom types
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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
