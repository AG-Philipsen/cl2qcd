/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "host_random.h"

#include <cstdio>
#include "../logger.hpp"
#include "../executables/exceptions.h"

#ifdef USE_PRNG_RANLUX
extern "C" {
#include "../ranlux/ranlxd.h"
}
#endif // USE_PRNG_RANLUX

/** Seed for the singleton random number generator rnd */
const unsigned long long int seed = 500000;

void prng_init(uint32_t seed)
{
#ifdef USE_PRNG_RANLUX
	// use maximum luxury level, should not be performance critical anyways
	rlxd_init(2, seed);
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
}

double prng_double()
{
#ifdef USE_PRNG_RANLUX
	double tmp;
	ranlxd(&tmp, 1);
	return tmp;
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
}

int prng_size()
{
#ifdef USE_PRNG_RANLUX
	return rlxd_size();
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
}

void prng_get(int* buf)
{
#ifdef USE_PRNG_RANLUX
	rlxd_get(buf);
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
}

void prng_set(int* buf)
{
#ifdef USE_PRNG_RANLUX
	rlxd_reset(buf);
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
}

void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha)
{
	hmc_float delta;
	hmc_float a0 ;
	hmc_float eta ;
#ifdef USE_PRNG_RANLUX
	do {
		double tmp[4];
		ranlxd(tmp, 4);
		delta = -log(tmp[0]) / alpha * pow(cos(2. * PI * tmp[1]), 2.) - log(tmp[2]) / alpha;
		a0 = 1. - delta;
		eta = tmp[3];
	} while ( (1. - 0.5 * delta) < eta * eta);
	double tmp[2];
	ranlxd(tmp, 2);
	hmc_float phi = 2.*PI * tmp[0];
	hmc_float theta = asin(2.*tmp[1] - 1.);
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
	dst[0] = a0;
	dst[1] = sqrt(1. - a0 * a0) * cos(theta) * cos(phi);
	dst[2] = sqrt(1. - a0 * a0) * cos(theta) * sin(phi);
	dst[3] = sqrt(1. - a0 * a0) * sin(theta);
}

void gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma)
{

	// SL: this fills real and imaginary part of a vector of "length" complex numbers
	//     with components drawn with a Gaussian distribution and variance sigma
	for(int idx = 0; idx < length; idx++) {
		gaussianNormalPair(&vector[idx].re, &vector[idx].im);
		vector[idx].re *= sigma;
		vector[idx].im *= sigma;
	}
	return;
	// SL: not yet tested
}

void gaussianNormalPair(hmc_float * z1, hmc_float * z2)
{
#ifdef USE_PRNG_RANLUX
	double tmp[2];
	ranlxd(tmp, 2);
	hmc_float u1 = 1.0 - tmp[0];
	hmc_float u2 = 1.0 - tmp[1];
#else // USE_PRNG_XXX
#error No implemented PRNG chosen
#endif // USE_PRNG_XXX
	hmc_float p  = sqrt(-2 * log(u1));
	*z1 = p * cos(2 * PI * u2);
	*z2 = p * sin(2 * PI * u2);
	// SL: not yet tested
}
