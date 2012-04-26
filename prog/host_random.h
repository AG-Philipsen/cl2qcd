
/** @file
 * Random number generations
 */

#ifndef _RANDOMH_
#define _RANDOMH_

#include <cstdlib>
#include <cstdio>
#include <sstream>

#include "globaldefs.h"
#include "types.h"
#include "host_geometry.h"
#include "host_use_timer.h"

/**
 * Seed the host prng.
 *
 * @param seed Seed for the underlying PRNG. Be aware of restrictions!
 */
void prng_init(uint32_t seed);

/**
 * Get a double precision random number from the generator
 *
 * @return A double precision number in [0, 1)
 */
double prng_double();

/** Storage type for state of the device random number generator */
typedef cl_ulong4 nr3_state_dev;

/**
 * init the given array with seeds for parallel random number generation on an OpenCL device
 *
 * @param[out] hmc_rndarray The array to feed with seeds
 * @param[in] file The file containing the binary seed
 * @param[in,out] inittime Timer to add execution time to
 */
void nr3_init_seeds(nr3_state_dev * const hmc_rndarray, char const * const seedfile, int const num_rndstates);

/** Construct new SU2 matrix using improved alg by Kennedy Pendleton */
void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha);

/**
 * Fill the real and imaginary parts of complex numbers in an
 * array with components drawn from a Gaussian distribution.
 *
 * \param[out] vector An array of complex numbers to write to
 * \param[in] length The amount of complex numbers to be copied
 * \param[in] sigma The variance of the Gaussian distribution
 */
void gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma);
/**
 * Generate two independent normal standard real numbers
 *
 * \param[out] z1 A real number
 * \param[out] z2 A real number
 */
void gaussianNormalPair(hmc_float * z1, hmc_float * z2);

#endif /* _RANDOMH_ */
