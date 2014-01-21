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

/**
 * @file operations on algebraelements
 */

#ifdef ENABLE_PRINTF
void print_ae(ae in)
{
	printf("%f, %f, %f, %f, %f, %f, %f, %f\n", in.e0, in.e1, in.e2, in.e3, in.e4, in.e5, in.e6, in.e7);
}
#endif

hmc_float ae_squarenorm(ae in)
{
	hmc_float result =
	  (in).e0 * (in).e0 +
	  (in).e1 * (in).e1 +
	  (in).e2 * (in).e2 +
	  (in).e3 * (in).e3 +
	  (in).e4 * (in).e4 +
	  (in).e5 * (in).e5 +
	  (in).e6 * (in).e6 +
	  (in).e7 * (in).e7;
	return result;
}

ae acc_factor_times_algebraelement(ae in, hmc_float factor, ae force_in)
{
	ae tmp;
	tmp.e0 = in.e0 + factor * force_in.e0;
	tmp.e1 = in.e1 + factor * force_in.e1;
	tmp.e2 = in.e2 + factor * force_in.e2;
	tmp.e3 = in.e3 + factor * force_in.e3;
	tmp.e4 = in.e4 + factor * force_in.e4;
	tmp.e5 = in.e5 + factor * force_in.e5;
	tmp.e6 = in.e6 + factor * force_in.e6;
	tmp.e7 = in.e7 + factor * force_in.e7;
	return tmp;
}

ae acc_algebraelement(ae in, ae force_in)
{
	ae tmp;
	tmp.e0 = in.e0 + force_in.e0;
	tmp.e1 = in.e1 + force_in.e1;
	tmp.e2 = in.e2 + force_in.e2;
	tmp.e3 = in.e3 + force_in.e3;
	tmp.e4 = in.e4 + force_in.e4;
	tmp.e5 = in.e5 + force_in.e5;
	tmp.e6 = in.e6 + force_in.e6;
	tmp.e7 = in.e7 + force_in.e7;
	return tmp;
}

// This function returns Tr[i * T_k * (in - in^\dag)], i.e. calculates the trace of
// i times a generator of the algebra su(3) times a 3x3-matrix that is (in - in^\dag).
// The result is stored in a su3-algebraelement:
//  - ae.e0 = Tr[i * T_0 * (in - in^\dag)]
//  - ae.e1 = Tr[i * T_1 * (in - in^\dag)]
//  - ae.e2 = Tr[i * T_2 * (in - in^\dag)]
// and so on (remember that T_k is the Gell Mann matrix lambda_k divided by 2).
ae tr_lambda_u(Matrix3x3 in)
{
	ae tmp;
	tmp.e0 = ( -in.e10.im - in.e01.im);
	tmp.e1 = (+in.e10.re - in.e01.re);
	tmp.e2 = (-in.e00.im + in.e11.im);
	tmp.e3 = (-in.e20.im - in.e02.im);
	tmp.e4 = (+in.e20.re - in.e02.re);
	tmp.e5 = (-in.e21.im - in.e12.im);
	tmp.e6 = (+in.e21.re - in.e12.re);
	tmp.e7 = (-in.e00.im - in.e11.im + 2.0 * in.e22.im) * 0.577350269189625;
	return tmp;
}

ae ae_times_factor(ae in, hmc_float factor)
{
	ae tmp;

	tmp.e0 = in.e0 * factor;
	tmp.e1 = in.e1 * factor;
	tmp.e2 = in.e2 * factor;
	tmp.e3 = in.e3 * factor;
	tmp.e4 = in.e4 * factor;
	tmp.e5 = in.e5 * factor;
	tmp.e6 = in.e6 * factor;
	tmp.e7 = in.e7 * factor;

	return tmp;
}

#ifdef GAUGEMOMENTASIZE_MEM // we only need the ifdef so we can use the non-field functions without having the field defined
ae getAe(__global const aeStorageType * const restrict in, const size_t idx)
{
#ifdef _USE_SOA_
	return (ae) {
		in[0 * GAUGEMOMENTA_STRIDE + idx],
		   in[1 * GAUGEMOMENTA_STRIDE + idx],
		   in[2 * GAUGEMOMENTA_STRIDE + idx],
		   in[3 * GAUGEMOMENTA_STRIDE + idx],
		   in[4 * GAUGEMOMENTA_STRIDE + idx],
		   in[5 * GAUGEMOMENTA_STRIDE + idx],
		   in[6 * GAUGEMOMENTA_STRIDE + idx],
		   in[7 * GAUGEMOMENTA_STRIDE + idx]
	};
#else
	return in[idx];
#endif
}

void putAe(__global aeStorageType * const restrict out, const size_t idx, const ae val)
{
#ifdef _USE_SOA_
	out[0 * GAUGEMOMENTA_STRIDE + idx] = val.e0;
	out[1 * GAUGEMOMENTA_STRIDE + idx] = val.e1;
	out[2 * GAUGEMOMENTA_STRIDE + idx] = val.e2;
	out[3 * GAUGEMOMENTA_STRIDE + idx] = val.e3;
	out[4 * GAUGEMOMENTA_STRIDE + idx] = val.e4;
	out[5 * GAUGEMOMENTA_STRIDE + idx] = val.e5;
	out[6 * GAUGEMOMENTA_STRIDE + idx] = val.e6;
	out[7 * GAUGEMOMENTA_STRIDE + idx] = val.e7;
#else
	out[idx] = val;
#endif
}

void update_gaugemomentum(const ae in, const hmc_float factor, const int global_link_pos, __global aeStorageType * const restrict out)
{
	ae tmp = getAe(out, global_link_pos);
	tmp = acc_factor_times_algebraelement(tmp, factor, in);
	putAe(out, global_link_pos, tmp);
}

#endif
