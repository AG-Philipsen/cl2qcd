/** @file
 * Operators und utility functions for the custom types
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "type_ops.hpp"

template<> void fill(hmc_complex* array, size_t num_elems)
{
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {static_cast<hmc_float>(i), static_cast<hmc_float>(-i)};
	}
}

template<> void fill(Matrixsu3* array, size_t num_elems)
{
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
			{ static_cast<hmc_float>(i + 1), static_cast<hmc_float>(i - 1) },
			{ static_cast<hmc_float>(i + 2), static_cast<hmc_float>(i - 2) },
			{ static_cast<hmc_float>(i + 3), static_cast<hmc_float>(i - 3) },
			{ static_cast<hmc_float>(i + 4), static_cast<hmc_float>(i - 4) },
			{ static_cast<hmc_float>(i + 5), static_cast<hmc_float>(i - 5) },
			{ static_cast<hmc_float>(i + 6), static_cast<hmc_float>(i - 6) },
			{ static_cast<hmc_float>(i + 7), static_cast<hmc_float>(i - 7) },
			{ static_cast<hmc_float>(i + 8), static_cast<hmc_float>(i - 8) },
			{ static_cast<hmc_float>(i + 9), static_cast<hmc_float>(i - 9) }
		};
	}
}
