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

/** @file
 * Device code implementing 3x3 matrices
 */

//operations_matrix.cl


#ifdef ENABLE_PRINTF
void print_matrix3x3(Matrix3x3 in)
{
	printf("(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n",
	       in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im,
	       in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im,
	       in.e20.re, in.e20.im, in.e21.re, in.e21.im, in.e22.re, in.e22.im);
	printf("\n");
}
#endif

inline Matrix3x3 zero_matrix3x3 ()
{
	return (Matrix3x3) {
		{0., 0.}, {0., 0.}, {0., 0.},
		{0., 0.}, {0., 0.}, {0., 0.},
		{0., 0.}, {0., 0.}, {0., 0.}
	};
}

inline Matrix3x3 identity_matrix3x3 ()
{
	return (Matrix3x3) {
		{1., 0.}, {0., 0.}, {0., 0.},
		{0., 0.}, {1., 0.}, {0., 0.},
		{0., 0.}, {0., 0.}, {1., 0.}
	};
}

inline Matrix3x3 multiply_matrix3x3_by_real (const Matrix3x3 in, const hmc_float factor)
{
	Matrix3x3 out = in;
	out.e00.re *= factor;
	out.e00.im *= factor;
	out.e01.re *= factor;
	out.e01.im *= factor;
	out.e02.re *= factor;
	out.e02.im *= factor;
	out.e10.re *= factor;
	out.e10.im *= factor;
	out.e11.re *= factor;
	out.e11.im *= factor;
	out.e12.re *= factor;
	out.e12.im *= factor;
	out.e20.re *= factor;
	out.e20.im *= factor;
	out.e21.re *= factor;
	out.e21.im *= factor;
	out.e22.re *= factor;
	out.e22.im *= factor;

	return out;
}

inline void scale_matrix3x3_by_real (Matrix3x3 * const restrict inout, const hmc_float factor)
{
	inout->e00.re *= factor;
	inout->e00.im *= factor;
	inout->e01.re *= factor;
	inout->e01.im *= factor;
	inout->e02.re *= factor;
	inout->e02.im *= factor;
	inout->e10.re *= factor;
	inout->e10.im *= factor;
	inout->e11.re *= factor;
	inout->e11.im *= factor;
	inout->e12.re *= factor;
	inout->e12.im *= factor;
	inout->e20.re *= factor;
	inout->e20.im *= factor;
	inout->e21.re *= factor;
	inout->e21.im *= factor;
	inout->e22.re *= factor;
	inout->e22.im *= factor;
}

inline Matrix3x3 multiply_matrix3x3_by_complex (Matrix3x3 in, hmc_complex factor)
{
	Matrix3x3 out;

	out.e00 = complexmult(in.e00, factor);
	out.e01 = complexmult(in.e01, factor);
	out.e02 = complexmult(in.e02, factor);
	out.e10 = complexmult(in.e10, factor);
	out.e11 = complexmult(in.e11, factor);
	out.e12 = complexmult(in.e12, factor);
	out.e20 = complexmult(in.e20, factor);
	out.e21 = complexmult(in.e21, factor);
	out.e22 = complexmult(in.e22, factor);

	return out;
}

inline Matrix3x3 multiply_matrix3x3 (const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e10.re + p.e02.re * q.e20.re
	             - p.e00.im * q.e00.im - p.e01.im * q.e10.im - p.e02.im * q.e20.im;
	out.e00.im = p.e00.re * q.e00.im + p.e01.re * q.e10.im + p.e02.re * q.e20.im
	             + p.e00.im * q.e00.re + p.e01.im * q.e10.re + p.e02.im * q.e20.re;

	out.e01.re = p.e00.re * q.e01.re + p.e01.re * q.e11.re + p.e02.re * q.e21.re
	             - p.e00.im * q.e01.im - p.e01.im * q.e11.im - p.e02.im * q.e21.im;
	out.e01.im = p.e00.re * q.e01.im + p.e01.re * q.e11.im + p.e02.re * q.e21.im
	             + p.e00.im * q.e01.re + p.e01.im * q.e11.re + p.e02.im * q.e21.re;

	out.e02.re = p.e00.re * q.e02.re + p.e01.re * q.e12.re + p.e02.re * q.e22.re
	             - p.e00.im * q.e02.im - p.e01.im * q.e12.im - p.e02.im * q.e22.im;
	out.e02.im = p.e00.re * q.e02.im + p.e01.re * q.e12.im + p.e02.re * q.e22.im
	             + p.e00.im * q.e02.re + p.e01.im * q.e12.re + p.e02.im * q.e22.re;

	out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e10.re + p.e12.re * q.e20.re
	             - p.e10.im * q.e00.im - p.e11.im * q.e10.im - p.e12.im * q.e20.im;
	out.e10.im = p.e10.re * q.e00.im + p.e11.re * q.e10.im + p.e12.re * q.e20.im
	             + p.e10.im * q.e00.re + p.e11.im * q.e10.re + p.e12.im * q.e20.re;

	out.e11.re = p.e10.re * q.e01.re + p.e11.re * q.e11.re + p.e12.re * q.e21.re
	             - p.e10.im * q.e01.im - p.e11.im * q.e11.im - p.e12.im * q.e21.im;
	out.e11.im = p.e10.re * q.e01.im + p.e11.re * q.e11.im + p.e12.re * q.e21.im
	             + p.e10.im * q.e01.re + p.e11.im * q.e11.re + p.e12.im * q.e21.re;

	out.e12.re = p.e10.re * q.e02.re + p.e11.re * q.e12.re + p.e12.re * q.e22.re
	             - p.e10.im * q.e02.im - p.e11.im * q.e12.im - p.e12.im * q.e22.im;
	out.e12.im = p.e10.re * q.e02.im + p.e11.re * q.e12.im + p.e12.re * q.e22.im
	             + p.e10.im * q.e02.re + p.e11.im * q.e12.re + p.e12.im * q.e22.re;

	out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e10.re + p.e22.re * q.e20.re
	             - p.e20.im * q.e00.im - p.e21.im * q.e10.im - p.e22.im * q.e20.im;
	out.e20.im = p.e20.re * q.e00.im + p.e21.re * q.e10.im + p.e22.re * q.e20.im
	             + p.e20.im * q.e00.re + p.e21.im * q.e10.re + p.e22.im * q.e20.re;

	out.e21.re = p.e20.re * q.e01.re + p.e21.re * q.e11.re + p.e22.re * q.e21.re
	             - p.e20.im * q.e01.im - p.e21.im * q.e11.im - p.e22.im * q.e21.im;
	out.e21.im = p.e20.re * q.e01.im + p.e21.re * q.e11.im + p.e22.re * q.e21.im
	             + p.e20.im * q.e01.re + p.e21.im * q.e11.re + p.e22.im * q.e21.re;

	out.e22.re = p.e20.re * q.e02.re + p.e21.re * q.e12.re + p.e22.re * q.e22.re
	             - p.e20.im * q.e02.im - p.e21.im * q.e12.im - p.e22.im * q.e22.im;
	out.e22.im = p.e20.re * q.e02.im + p.e21.re * q.e12.im + p.e22.re * q.e22.im
	             + p.e20.im * q.e02.re + p.e21.im * q.e12.re + p.e22.im * q.e22.re;

	return out;
}

inline Matrix3x3 multiply_matrix3x3_dagger(const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;

	out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
	             + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
	out.e00.im = -p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
	             + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

	out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
	             + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
	out.e01.im = -p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
	             + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;

	out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
	             + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
	out.e02.im = -p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
	             + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;

	out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
	             + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
	out.e10.im = -p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
	             + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;

	out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
	             + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
	out.e11.im = -p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
	             + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;

	out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
	             + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
	out.e12.im = -p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
	             + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;

	out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
	             + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
	out.e20.im = -p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
	             + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;

	out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
	             + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
	out.e21.im = -p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
	             + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;

	out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
	             + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
	out.e22.im = -p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
	             + p.e20.im * q.e20.re + p.e21.im * q.e21.re + p.e22.im * q.e22.re;
	return out;
}

//this returns p^dagger*q^dagger
inline Matrix3x3 multiply_matrix3x3_dagger_dagger(const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;

	out.e00.re = p.e00.re * q.e00.re + p.e10.re * q.e01.re + p.e20.re * q.e02.re
	             - p.e00.im * q.e00.im - p.e10.im * q.e01.im - p.e20.im * q.e02.im;
	out.e00.im = -p.e00.re * q.e00.im - p.e10.re * q.e01.im - p.e20.re * q.e02.im
	             - p.e00.im * q.e00.re - p.e10.im * q.e01.re - p.e20.im * q.e02.re;

	out.e01.re = p.e00.re * q.e10.re + p.e10.re * q.e11.re + p.e20.re * q.e12.re
	             - p.e00.im * q.e10.im - p.e10.im * q.e11.im - p.e20.im * q.e12.im;
	out.e01.im = -p.e00.re * q.e10.im - p.e10.re * q.e11.im - p.e20.re * q.e12.im
	             - p.e00.im * q.e10.re - p.e10.im * q.e11.re - p.e20.im * q.e12.re;

	out.e02.re = p.e00.re * q.e20.re + p.e10.re * q.e21.re + p.e20.re * q.e22.re
	             - p.e00.im * q.e20.im - p.e10.im * q.e21.im - p.e20.im * q.e22.im;
	out.e02.im = -p.e00.re * q.e20.im - p.e10.re * q.e21.im - p.e20.re * q.e22.im
	             - p.e00.im * q.e20.re - p.e10.im * q.e21.re - p.e20.im * q.e22.re;

	out.e10.re = p.e01.re * q.e00.re + p.e11.re * q.e01.re + p.e21.re * q.e02.re
	             - p.e01.im * q.e00.im - p.e11.im * q.e01.im - p.e21.im * q.e02.im;
	out.e10.im = -p.e01.re * q.e00.im - p.e11.re * q.e01.im - p.e21.re * q.e02.im
	             - p.e01.im * q.e00.re - p.e11.im * q.e01.re - p.e21.im * q.e02.re;

	out.e11.re = p.e01.re * q.e10.re + p.e11.re * q.e11.re + p.e21.re * q.e12.re
	             - p.e01.im * q.e10.im - p.e11.im * q.e11.im - p.e21.im * q.e12.im;
	out.e11.im = -p.e01.re * q.e10.im - p.e11.re * q.e11.im - p.e21.re * q.e12.im
	             - p.e01.im * q.e10.re - p.e11.im * q.e11.re - p.e21.im * q.e12.re;

	out.e12.re = p.e01.re * q.e20.re + p.e11.re * q.e21.re + p.e21.re * q.e22.re
	             - p.e01.im * q.e20.im - p.e11.im * q.e21.im - p.e21.im * q.e22.im;
	out.e12.im = -p.e01.re * q.e20.im - p.e11.re * q.e21.im - p.e21.re * q.e22.im
	             - p.e01.im * q.e20.re - p.e11.im * q.e21.re - p.e21.im * q.e22.re;

	out.e20.re = p.e02.re * q.e00.re + p.e12.re * q.e01.re + p.e22.re * q.e02.re
	             - p.e02.im * q.e00.im - p.e12.im * q.e01.im - p.e22.im * q.e02.im;
	out.e20.im = -p.e02.re * q.e00.im - p.e12.re * q.e01.im - p.e22.re * q.e02.im
	             - p.e02.im * q.e00.re - p.e12.im * q.e01.re - p.e22.im * q.e02.re;

	out.e21.re = p.e02.re * q.e10.re + p.e12.re * q.e11.re + p.e22.re * q.e12.re
	             - p.e02.im * q.e10.im - p.e12.im * q.e11.im - p.e22.im * q.e12.im;
	out.e21.im = -p.e02.re * q.e10.im - p.e12.re * q.e11.im - p.e22.re * q.e12.im
	             - p.e02.im * q.e10.re - p.e12.im * q.e11.re - p.e22.im * q.e12.re;

	out.e22.re = p.e02.re * q.e20.re + p.e12.re * q.e21.re + p.e22.re * q.e22.re
	             - p.e02.im * q.e20.im - p.e12.im * q.e21.im - p.e22.im * q.e22.im;
	out.e22.im = -p.e02.re * q.e20.im - p.e12.re * q.e21.im - p.e22.re * q.e22.im
	             - p.e02.im * q.e20.re - p.e12.im * q.e21.re - p.e22.im * q.e22.re;





	/*


	  out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
	               + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
	  out.e00.im = -p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
	               + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

	  out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
	               + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
	  out.e01.im = -p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
	               + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;

	  out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
	               + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
	  out.e02.im = -p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
	               + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;

	  out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
	               + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
	  out.e10.im = -p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
	               + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;

	  out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
	               + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
	  out.e11.im = -p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
	               + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;

	  out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
	               + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
	  out.e12.im = -p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
	               + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;

	  out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
	               + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
	  out.e20.im = -p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
	               + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;

	  out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
	               + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
	  out.e21.im = -p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
	               + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;

	  out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
	               + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
	  out.e22.im = -p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
	               + p.e20.im * q.e20.re + p.e21.im * q.e21.re + p.e22.im * q.e22.re;
	*/
	return out;
}


inline Matrix3x3 add_matrix3x3 ( const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re + q.e00.re;
	out.e00.im = p.e00.im + q.e00.im;
	out.e01.re = p.e01.re + q.e01.re;
	out.e01.im = p.e01.im + q.e01.im;
	out.e02.re = p.e02.re + q.e02.re;
	out.e02.im = p.e02.im + q.e02.im;

	out.e10.re = p.e10.re + q.e10.re;
	out.e10.im = p.e10.im + q.e10.im;
	out.e11.re = p.e11.re + q.e11.re;
	out.e11.im = p.e11.im + q.e11.im;
	out.e12.re = p.e12.re + q.e12.re;
	out.e12.im = p.e12.im + q.e12.im;

	out.e20.re = p.e20.re + q.e20.re;
	out.e20.im = p.e20.im + q.e20.im;
	out.e21.re = p.e21.re + q.e21.re;
	out.e21.im = p.e21.im + q.e21.im;
	out.e22.re = p.e22.re + q.e22.re;
	out.e22.im = p.e22.im + q.e22.im;

	return out;
}

inline void aggregate_matrix3x3(Matrix3x3 * const restrict p, const Matrix3x3 q)
{
	p->e00.re += q.e00.re;
	p->e00.im += q.e00.im;
	p->e01.re += q.e01.re;
	p->e01.im += q.e01.im;
	p->e02.re += q.e02.re;
	p->e02.im += q.e02.im;

	p->e10.re += q.e10.re;
	p->e10.im += q.e10.im;
	p->e11.re += q.e11.re;
	p->e11.im += q.e11.im;
	p->e12.re += q.e12.re;
	p->e12.im += q.e12.im;

	p->e20.re += q.e20.re;
	p->e20.im += q.e20.im;
	p->e21.re += q.e21.re;
	p->e21.im += q.e21.im;
	p->e22.re += q.e22.re;
	p->e22.im += q.e22.im;
}

inline Matrix3x3 subtract_matrix3x3 ( const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re - q.e00.re;
	out.e00.im = p.e00.im - q.e00.im;
	out.e01.re = p.e01.re - q.e01.re;
	out.e01.im = p.e01.im - q.e01.im;
	out.e02.re = p.e02.re - q.e02.re;
	out.e02.im = p.e02.im - q.e02.im;

	out.e10.re = p.e10.re - q.e10.re;
	out.e10.im = p.e10.im - q.e10.im;
	out.e11.re = p.e11.re - q.e11.re;
	out.e11.im = p.e11.im - q.e11.im;
	out.e12.re = p.e12.re - q.e12.re;
	out.e12.im = p.e12.im - q.e12.im;

	out.e20.re = p.e20.re - q.e20.re;
	out.e20.im = p.e20.im - q.e20.im;
	out.e21.re = p.e21.re - q.e21.re;
	out.e21.im = p.e21.im - q.e21.im;
	out.e22.re = p.e22.re - q.e22.re;
	out.e22.im = p.e22.im - q.e22.im;

	return out;
}

inline Matrix3x3 subtract_matrix3x3_dagger (const Matrix3x3 p, const Matrix3x3 q)
{
	Matrix3x3 out;

	out.e00.re = p.e00.re - q.e00.re;
	out.e00.im = p.e00.im + q.e00.im;
	out.e01.re = p.e01.re - q.e10.re;
	out.e01.im = p.e01.im + q.e10.im;
	out.e02.re = p.e02.re - q.e20.re;
	out.e02.im = p.e02.im + q.e20.im;

	out.e10.re = p.e10.re - q.e01.re;
	out.e10.im = p.e10.im + q.e01.im;
	out.e11.re = p.e11.re - q.e11.re;
	out.e11.im = p.e11.im + q.e11.im;
	out.e12.re = p.e12.re - q.e21.re;
	out.e12.im = p.e12.im + q.e21.im;

	out.e20.re = p.e20.re - q.e02.re;
	out.e20.im = p.e20.im + q.e02.im;
	out.e21.re = p.e21.re - q.e12.re;
	out.e21.im = p.e21.im + q.e12.im;
	out.e22.re = p.e22.re - q.e22.re;
	out.e22.im = p.e22.im + q.e22.im;

	return out;
}

inline Matrix3x3 copy_matrix3x3(const Matrix3x3 p)
{
	Matrix3x3 out;
	out.e00.re = p.e00.re;
	out.e00.im = p.e00.im;
	out.e01.re = p.e01.re ;
	out.e01.im = p.e01.im;
	out.e02.re = p.e02.re;
	out.e02.im = p.e02.im;

	out.e10.re = p.e10.re;
	out.e10.im = p.e10.im;
	out.e11.re = p.e11.re;
	out.e11.im = p.e11.im;
	out.e12.re = p.e12.re;
	out.e12.im = p.e12.im;

	out.e20.re = p.e20.re;
	out.e20.im = p.e20.im;
	out.e21.re = p.e21.re;
	out.e21.im = p.e21.im;
	out.e22.re = p.e22.re;
	out.e22.im = p.e22.im;

	return out;
}

inline hmc_complex trace_matrix3x3 (const Matrix3x3 p)
{
	hmc_complex out;
	out.re = p.e00.re;
	out.im = p.e00.im;
	out.re += p.e11.re;
	out.im += p.e11.im;
	out.re += p.e22.re;
	out.im += p.e22.im;

	return out;
}

inline Matrix3x3 adjoint_matrix3x3 (const Matrix3x3 p)
{

	Matrix3x3 out;
	out.e00.re = p.e00.re;
	out.e00.im = - p.e00.im;
	out.e01.re = p.e10.re;
	out.e01.im = - p.e10.im;
	out.e02.re = p.e20.re;
	out.e02.im = - p.e20.im;


	out.e10.re = p.e01.re;
	out.e10.im = - p.e01.im;
	out.e11.re = p.e11.re;
	out.e11.im = - p.e11.im;
	out.e12.re = p.e21.re;
	out.e12.im = - p.e21.im;

	out.e20.re = p.e02.re;
	out.e20.im = - p.e02.im;
	out.e21.re = p.e12.re;
	out.e21.im = - p.e12.im;
	out.e22.re = p.e22.re;
	out.e22.im = - p.e22.im;

	return out;

}

inline Matrix3x3 matrix_su3to3x3 (const Matrixsu3 p)
{
	return (Matrix3x3) {
		p.e00, p.e01, p.e02,
		       p.e10, p.e11, p.e12,
		       p.e20, p.e21, p.e22
	};
}


inline Matrixsu3 matrix_3x3tosu3 (const Matrix3x3 p)
{
	return (Matrixsu3) {
		p.e00, p.e01, p.e02,
		       p.e10, p.e11, p.e12,
		       p.e20, p.e21, p.e22
	};
}

inline hmc_float absoluteDifference_matrix3x3(Matrix3x3 mat1, Matrix3x3 mat2)
{
	hmc_float result = 0.0;

	result += fabs((mat1).e00.re - (mat2).e00.re);
	result += fabs((mat1).e00.im - (mat2).e00.im);
	result += fabs((mat1).e01.re - (mat2).e01.re);
	result += fabs((mat1).e01.im - (mat2).e01.im);
	result += fabs((mat1).e02.re - (mat2).e02.re);
	result += fabs((mat1).e02.im - (mat2).e02.im);
	result += fabs((mat1).e10.re - (mat2).e10.re);
	result += fabs((mat1).e10.im - (mat2).e10.im);
	result += fabs((mat1).e11.re - (mat2).e11.re);
	result += fabs((mat1).e11.im - (mat2).e11.im);
	result += fabs((mat1).e12.re - (mat2).e12.re);
	result += fabs((mat1).e12.im - (mat2).e12.im);
	result += fabs((mat1).e20.re - (mat2).e20.re);
	result += fabs((mat1).e20.im - (mat2).e20.im);
	result += fabs((mat1).e21.re - (mat2).e21.re);
	result += fabs((mat1).e21.im - (mat2).e21.im);
	result += fabs((mat1).e22.re - (mat2).e22.re);
	result += fabs((mat1).e22.im - (mat2).e22.im);

	return result;
}

inline Matrixsu3 project_anti_herm(Matrix3x3 in)
{
	Matrixsu3 tmp;
	hmc_float tr_omega ;

	//in is not a su3 matrix, so reconstruct does not have to be applied here!!
	tr_omega = (in.e00.im + in.e11.im + in.e22.im) / 3.0 ;

	tmp.e00.re = 0.0 ;
	tmp.e00.im = -in.e00.im + tr_omega ;

	tmp.e11.re = 0.0 ;
	tmp.e11.im = -in.e11.im + tr_omega ;

	tmp.e22.re = 0.0 ;
	tmp.e22.im = -in.e22.im + tr_omega ;

	tmp.e01.re = in.e01.re - in.e10.re ;
	tmp.e01.re /= -2.0 ;
	tmp.e01.im = in.e01.im + in.e10.im ;
	tmp.e01.im /= -2.0 ;
	tmp.e10.re  = -tmp.e01.re ;
	tmp.e10.im  =  tmp.e01.im ;

	tmp.e02.re = in.e02.re - in.e20.re ;
	tmp.e02.re /= -2.0 ;
	tmp.e02.im = in.e02.im + in.e20.im ;
	tmp.e02.im /= -2.0 ;
	tmp.e20.re  = -tmp.e02.re ;
	tmp.e20.im  =  tmp.e02.im ;
	tmp.e12.re = in.e12.re - in.e21.re ;
	tmp.e12.re /= -2.0 ;
	tmp.e12.im = in.e12.im + in.e21.im ;
	tmp.e12.im /= -2.0 ;
	tmp.e21.re  = -tmp.e12.re ;
	tmp.e21.im  =  tmp.e12.im ;
	// flip sign, this is just taken from the tmlqcd analogue
	tmp = multiply_matrixsu3_by_real(tmp, -1.0) ;

	return tmp;

}

//This function returns the traceless antihermitian part of the input Matrix3x3.
//It is nothing but
//  M_TA = 1/2(M-M^dag) - 1/6 Tr(M-M^dag) Id
//where Id is the 3x3 identity matrix.
//NOTE: This function is exactely the same as project_anti_herm but
//      here the sign flip in the end doesn't appear because it is
//      done during the function (the factor of the division is +2.0 and not -2.0).
//      I checked in the tmLQCD / smearing / utils_project_antiherm.c file
//      and the function void project_antiherm(su3 *omega) is exactly this.
//NOTE: The output is an ANTIhermitean matrix so it is NOT a Matrixsu3 but only a Matrix3x3.
/**
 * @todo Fix this doubling of the code, basically figuring out whether it is correct
 *       (I do not think so) that project_anti_herm returns a Matrixsu3 object. In case
 *       substitute it by traceless_antihermitian_part.
 */
inline Matrix3x3 traceless_antihermitian_part(Matrix3x3 in)
{
	Matrix3x3 out;
	hmc_float trace = (in.e00.im + in.e11.im + in.e22.im) / 3.0;
	
	out.e00.re = 0.0;
	out.e11.re = 0.0;
	out.e22.re = 0.0;
	out.e00.im = in.e00.im - trace;
	out.e11.im = in.e11.im - trace;
	out.e22.im = in.e22.im - trace;
	
	out.e01.re = 0.5 * (in.e01.re - in.e10.re);
	out.e01.im = 0.5 * (in.e01.im + in.e10.im);
	out.e02.re = 0.5 * (in.e02.re - in.e20.re);
	out.e02.im = 0.5 * (in.e02.im + in.e20.im);
	out.e12.re = 0.5 * (in.e12.re - in.e21.re);
	out.e12.im = 0.5 * (in.e12.im + in.e21.im);
	
	out.e10.re = -out.e01.re;
	out.e10.im =  out.e01.im;
	out.e20.re = -out.e02.re;
	out.e20.im =  out.e02.im;
	out.e21.re = -out.e12.re;
	out.e21.im =  out.e12.im;
	
	return out;
}
