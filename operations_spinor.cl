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
 * Device code for operations on spinors
 */

#ifdef ENABLE_PRINTF
void print_spinor(spinor in)
{
	print_su3vec(in.e0);
	printf("\n");
	print_su3vec(in.e1);
	printf("\n");
	print_su3vec(in.e2);
	printf("\n");
	print_su3vec(in.e3);
	printf("\n");
}
#endif

spinor set_spinor_zero()
{
	spinor tmp;
	tmp.e0 = set_su3vec_zero();
	tmp.e1 = set_su3vec_zero();
	tmp.e2 = set_su3vec_zero();
	tmp.e3 = set_su3vec_zero();
	return tmp;
}

spinor set_spinor_cold()
{
	spinor tmp;
	tmp.e0 = set_su3vec_cold();
	tmp.e1 = set_su3vec_cold();
	tmp.e2 = set_su3vec_cold();
	tmp.e3 = set_su3vec_cold();
	return tmp;
}

hmc_float spinor_squarenorm(spinor in)
{
	hmc_float res = 0;
	res += su3vec_squarenorm(in.e0);
	res += su3vec_squarenorm(in.e1);
	res += su3vec_squarenorm(in.e2);
	res += su3vec_squarenorm(in.e3);
	return res;
}

hmc_complex spinor_scalarproduct(spinor in1, spinor in2)
{
	hmc_complex res = hmc_complex_zero;
	hmc_complex tmp;
	tmp = su3vec_scalarproduct(in1.e0, in2.e0);
	res.re += tmp.re;
	res.im += tmp.im;
	tmp = su3vec_scalarproduct(in1.e1, in2.e1);
	res.re += tmp.re;
	res.im += tmp.im;
	tmp = su3vec_scalarproduct(in1.e2, in2.e2);
	res.re += tmp.re;
	res.im += tmp.im;
	tmp = su3vec_scalarproduct(in1.e3, in2.e3);
	res.re += tmp.re;
	res.im += tmp.im;
	return res;
}

spinor real_multiply_spinor(spinor in, hmc_float factor)
{
	spinor tmp;
	tmp.e0 = su3vec_times_real(in.e0, factor);
	tmp.e1 = su3vec_times_real(in.e1, factor);
	tmp.e2 = su3vec_times_real(in.e2, factor);
	tmp.e3 = su3vec_times_real(in.e3, factor);
	return tmp;
}

spinor spinor_times_complex(spinor in, hmc_complex factor)
{
	spinor tmp;
	tmp.e0 = su3vec_times_complex(in.e0, factor);
	tmp.e1 = su3vec_times_complex(in.e1, factor);
	tmp.e2 = su3vec_times_complex(in.e2, factor);
	tmp.e3 = su3vec_times_complex(in.e3, factor);
	return tmp;
}

spinor spinor_times_complex_conj(spinor in, hmc_complex factor)
{
	spinor tmp;
	tmp.e0 = su3vec_times_complex_conj(in.e0, factor);
	tmp.e1 = su3vec_times_complex_conj(in.e1, factor);
	tmp.e2 = su3vec_times_complex_conj(in.e2, factor);
	tmp.e3 = su3vec_times_complex_conj(in.e3, factor);
	return tmp;
}

spinor spinor_dim(spinor in1, spinor in2)
{
	spinor tmp;
	tmp.e0 = su3vec_dim(in1.e0, in2.e0);
	tmp.e1 = su3vec_dim(in1.e1, in2.e1);
	tmp.e2 = su3vec_dim(in1.e2, in2.e2);
	tmp.e3 = su3vec_dim(in1.e3, in2.e3);
	return tmp;
}

spinor spinor_acc(spinor in1, spinor in2)
{
	spinor tmp;
	tmp.e0 = su3vec_acc(in1.e0, in2.e0);
	tmp.e1 = su3vec_acc(in1.e1, in2.e1);
	tmp.e2 = su3vec_acc(in1.e2, in2.e2);
	tmp.e3 = su3vec_acc(in1.e3, in2.e3);
	return tmp;
}

spinor spinor_acc_acc(spinor in1, spinor in2, spinor in3)
{
	spinor tmp;
	tmp.e0 = su3vec_acc_acc(in1.e0, in2.e0, in3.e0);
	tmp.e1 = su3vec_acc_acc(in1.e1, in2.e1, in3.e1);
	tmp.e2 = su3vec_acc_acc(in1.e2, in2.e2, in3.e2);
	tmp.e3 = su3vec_acc_acc(in1.e3, in2.e3, in3.e3);
	return tmp;
}

