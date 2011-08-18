/** @file
 * Random number generations
 */

#ifndef _RANDOMH_
#define _RANDOMH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_geometry.h"
#include "host_use_timer.h"

/** Utility 64-bit integer for random number generation */
typedef unsigned long long int Ullong;
/** Utility 32-bit integer for random number generation */
typedef unsigned int Uint;

/** Seed for the singleton random number generator rnd */
const unsigned long long int seed = 500000;

/** The Random number generator described in Numerical Recipes 3 */
struct Random {

	/** Random number state */
	Ullong u,v,w;

	/**
	 * Initializes the random number generator.
	 *
	 * @param j Seed for the random number generator state
	 */
	Random(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v;
		int64();
		v = u;
		int64();
		w = v;
		int64();
	}

	/**
	 * Generate a random 64-bit integer.
	 */
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17;
		v ^= v << 31;
		v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21);
		x ^= x >> 35;
		x ^= x << 4;
		return (x + v) ^ w;
	}
	/**
	 * Generate a random 64-bit floating point number.
	 */
	inline double doub() {
		return 5.42101086242752217E-20 * int64();
	}
	/**
	 * Generate a random 32-bit integer.
	 */
	inline Uint int32() {
		return (Uint)int64();
	}
};

/** The singleton single-threaded random number generator */
extern Random rnd;

/**
 * Get 1,2,3 in random order
 *
 * @param[out] rand Storage location for the result
 */
void random_1_2_3 (int rand[3]);

/**
 * init the given array with seeds for parallel random number generation on an OpenCL device
 *
 * @param[out] hmc_rndarray The array to feed with seeds
 * @param[in] file The file containing the binary seed
 * @param[in,out] inittime Timer to add execution time to
 * @return Error code as defined in hmcerrs.h
 *         @li HMC_FILEERROR    if the seeding file cannot be opened
 *         @li HMC_INVALIDVALUE if the seeding file does not contain enough bytes
 *         @li HMC_SUCCESS      otherwise
 */
int init_random_seeds(hmc_ocl_ran * const hmc_rndarray, char const * const seedfile);

/** Construct new SU2 matrix using improved alg by Kennedy Pendleton */
void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha);

#endif /* _RANDOMH_ */
