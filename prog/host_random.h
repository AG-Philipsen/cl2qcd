
/** @file
 * Random number generations
 */

#ifndef _RANDOMH_
#define _RANDOMH_

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

int prng_size();

void prng_get(int* buf);

void prng_set(int* buf);

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
