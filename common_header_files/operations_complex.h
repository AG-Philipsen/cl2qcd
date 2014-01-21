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

#ifndef OPERATIONS_COMPLEX
#define OPERATIONS_COMPLEX

inline hmc_complex complexconj(hmc_complex in)
{
	in.im = -(in.im);
	return in;
}

inline hmc_complex complexmult(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re * b.re - a.im * b.im;
	res.im = a.im * b.re + a.re * b.im;
	return res;
}

inline hmc_complex complexadd(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re + b.re;
	res.im = a.im + b.im;
	return res;
}

inline hmc_complex complexsubtract(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re - b.re;
	res.im = a.im - b.im;
	return res;
}

inline hmc_complex complexdivide(const hmc_complex numerator, const hmc_complex denominator)
{
	hmc_float norm = denominator.re * denominator.re + denominator.im * denominator.im;
	hmc_complex res;
	res.re = (numerator.re * denominator.re + numerator.im * denominator.im ) / norm;
	res.im = (numerator.im * denominator.re - numerator.re * denominator.im ) / norm;
	return res;
}

#ifdef _INKERNEL_


/** This type can be used to store hmc_complex as it has the same size in bytes */
#ifdef _USEDOUBLEPREC_
typedef float4 hmc_complex_store_type;
#else
typedef float2 hmc_complex_store_type;
#endif

/**
 * Workaround for complex constatns not being loaded properly on APP 2.5
 *
 * For a further discussion of this bug see
 * http://forums.amd.com/devforum/messageview.cfm?catid=390&threadid=156057&enterthread=y
 *
 * @todo make this only be used on APP 2.5
 */
inline hmc_complex complexLoadHack(__global const hmc_complex * p)
{
#ifdef USE_CPLX_LOAD_HACK
	union {
		hmc_complex_store_type v;
		hmc_complex c;
	} tmp;
	tmp.v = *((__global const hmc_complex_store_type*) p);
	return tmp.c;
#else
	return p[0];
#endif
}

#endif // ifdef _INKERNEL_

#endif /* OPERATIONS_COMPLEX */
