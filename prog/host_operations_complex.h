/** @file
 * Mathematical operations on the hmc_complex type.
 *
 * @todo does it really make sense to use pointers for the args?
 * @todo it would probably make sense to allow inlining of these functions (high invocation density)
 */

#ifndef _OPERATIONS_COMPLEXH_
#define _OPERATIONS_COMPLEXH_

#include <iostream>
#include <math.h>
#include "host_random.h"
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"

/**
 * Calculate the complex conjugate.
 *
 * \param in A complex number
 * \return The complex conjugate of in
 */
hmc_complex complexconj(hmc_complex *in); 
/**
 * Multiply two complex numbers
 *
 * \param a A complex number
 * \param b A complex number
 * \return The result of multiplying a and b
 */
hmc_complex complexmult(hmc_complex *a, hmc_complex *b); 
/**
 * Add two complex numbers
 *
 * \param a A complex number
 * \param b A complex number
 * \return The result of adding a and b
 */
hmc_complex complexadd(hmc_complex *a, hmc_complex *b);
/**
 * Substract two complex numbers
 *
 * \param a A complex number
 * \param b A complex number
 * \return The result of subtracting b from a
 */
hmc_complex complexsubtract(hmc_complex *a, hmc_complex *b); 
/**
 * Accumulates complex numbers.
 *
 * Given an accumulation variable and an increment adds the
 * increment to the accumulation variable on each invocation.
 *
 * \param[in,out] inout A complex number used as an accumulator
 * \param[in] incr A complex number specifying the increment
 * \return Error state as defined in hmcerrs.h
 */
hmc_error complexaccumulate(hmc_complex *inout, hmc_complex *incr); 
/**
 * Divide two complex numbers
 *
 * \param numerator A complex number to be divided
 * \param denominator A complex number to divide by
 * \return The result of dividing numerator by denominator
 */
hmc_complex complexdivide(hmc_complex* numerator, hmc_complex* denominator);
/**
 * Copy an array of complex numbers.
 *
 * \param[in] source Pointer to the complex number to be copied
 * \param[out] dest Pointer to copy the complex numbers to
 * \param[in] length The amount of complex numbers to be copied
 * \return Error state as defined in hmcerrs.h
 */
hmc_error complexcopy(hmc_complex* source, hmc_complex* dest, int length);
/**
 * Scale a complex number in-place
 *
 * \param[in,out] a The complex number to scale
 * \param[in] b A real scaling factor (may be negative)
 * \return Error state as defined in hmcerrs.h
 */
hmc_error complexmult_real(hmc_complex *a, hmc_float *b); 
/**
 * Fill the real and imaginary parts of complex numbers in an
 * array with components drawn from a Gaussian distribution.
 *
 * \param[out] vector An array of complex numbers to write to
 * \param[in] length The amount of complex numbers to be copied
 * \param[in] sigma The variance of the Gaussian distribution
 * \return Error state as defined in hmcerrs.h
 */
hmc_error gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma);
/**
 * Generate two independent normal standard real numbers
 *
 * \param[out] z1 A real number
 * \param[out] z2 A real number
 */
void gaussianNormalPair(hmc_float * z1, hmc_float * z2);

#endif

