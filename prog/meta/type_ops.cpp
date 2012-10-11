/** @file
 * Operators und utility functions for the custom types
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "type_ops.hpp"

template<> void fill(hmc_complex* array, size_t num_elems, int seed)
{
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {static_cast<hmc_float>(i) + seed, static_cast<hmc_float>(-i) + seed};
	}
}

template<> void fill(Matrixsu3* array, size_t num_elems, int seed)
{
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
			{ static_cast<hmc_float>(i + 1) + seed, static_cast<hmc_float>(i - 1) + seed },
			{ static_cast<hmc_float>(i + 2) + seed, static_cast<hmc_float>(i - 2) + seed },
			{ static_cast<hmc_float>(i + 3) + seed, static_cast<hmc_float>(i - 3) + seed },
			{ static_cast<hmc_float>(i + 4) + seed, static_cast<hmc_float>(i - 4) + seed },
			{ static_cast<hmc_float>(i + 5) + seed, static_cast<hmc_float>(i - 5) + seed },
			{ static_cast<hmc_float>(i + 6) + seed, static_cast<hmc_float>(i - 6) + seed },
			{ static_cast<hmc_float>(i + 7) + seed, static_cast<hmc_float>(i - 7) + seed },
			{ static_cast<hmc_float>(i + 8) + seed, static_cast<hmc_float>(i - 8) + seed },
			{ static_cast<hmc_float>(i + 9) + seed, static_cast<hmc_float>(i - 9) + seed }
		};
	}
}

template<> void fill(spinor* array, size_t num_elems, int seed)
{
	for(size_t i = 0; i < num_elems; i++) {
		array[i] = {
			{
				{ static_cast<hmc_float>(i + 1) + seed, static_cast<hmc_float>(i - 1) + seed },
				{ static_cast<hmc_float>(i + 2) + seed, static_cast<hmc_float>(i - 2) + seed },
				{ static_cast<hmc_float>(i + 3) + seed, static_cast<hmc_float>(i - 3) + seed },
			}, {
				{ static_cast<hmc_float>(i + 4) + seed, static_cast<hmc_float>(i - 4) + seed },
				{ static_cast<hmc_float>(i + 5) + seed, static_cast<hmc_float>(i - 5) + seed },
				{ static_cast<hmc_float>(i + 6) + seed, static_cast<hmc_float>(i - 6) + seed },
			}, {
				{ static_cast<hmc_float>(i + 7) + seed, static_cast<hmc_float>(i - 7) + seed },
				{ static_cast<hmc_float>(i + 8) + seed, static_cast<hmc_float>(i - 8) + seed },
				{ static_cast<hmc_float>(i + 9) + seed, static_cast<hmc_float>(i - 9) + seed }
			}, {
				{ static_cast<hmc_float>(i + 10) + seed, static_cast<hmc_float>(i - 10) + seed },
				{ static_cast<hmc_float>(i + 11) + seed, static_cast<hmc_float>(i - 11) + seed },
				{ static_cast<hmc_float>(i + 12) + seed, static_cast<hmc_float>(i - 12) + seed }
			}
		};
	}
}
