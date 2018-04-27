/*
 * Copyright (c) 2011,2013 Matthias Bach
 * Copyright (c) 2011 Christopher Pinke
 * Copyright (c) 2013,2014,2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * Device code for operations on SU(3) vectors
 *
 */

#ifdef ENABLE_PRINTF
void print_su3vec(su3vec in)
{
    printf("(%f,%f)\t(%f,%f)\t(%f,%f)", in.e0.re, in.e0.im, in.e1.re, in.e1.im, in.e2.re, in.e2.im);
}
#endif

hmc_float su3vec_squarenorm(su3vec in)
{
    return in.e0.re * in.e0.re + in.e0.im * in.e0.im + in.e1.re * in.e1.re + in.e1.im * in.e1.im + in.e2.re * in.e2.re +
           in.e2.im * in.e2.im;
}

hmc_complex su3vec_scalarproduct(su3vec in1, su3vec in2)
{
    hmc_complex tmp;
    tmp.re = in1.e0.re * in2.e0.re + in1.e0.im * in2.e0.im + in1.e1.re * in2.e1.re + in1.e1.im * in2.e1.im +
             in1.e2.re * in2.e2.re + in1.e2.im * in2.e2.im;
    tmp.im = in1.e0.re * in2.e0.im - in2.e0.re * in1.e0.im + in1.e1.re * in2.e1.im - in2.e1.re * in1.e1.im +
             in1.e2.re * in2.e2.im - in2.e2.re * in1.e2.im;
    return tmp;
}

hmc_float su3vec_scalarproduct_real_part(su3vec in1, su3vec in2)
{
    hmc_float tmp;
    tmp = in1.e0.re * in2.e0.re + in1.e0.im * in2.e0.im + in1.e1.re * in2.e1.re + in1.e1.im * in2.e1.im +
          in1.e2.re * in2.e2.re + in1.e2.im * in2.e2.im;
    return tmp;
}

su3vec set_su3vec_zero()
{
    su3vec tmp;
    (tmp).e0 = hmc_complex_zero;
    (tmp).e1 = hmc_complex_zero;
    (tmp).e2 = hmc_complex_zero;
    return tmp;
}

su3vec set_su3vec_cold()
{
    su3vec tmp;
    (tmp).e0 = hmc_complex_one;
    (tmp).e1 = hmc_complex_one;
    (tmp).e2 = hmc_complex_one;
    return tmp;
}

su3vec su3vec_times_real(su3vec in, hmc_float factor)
{
    su3vec tmp;
    tmp.e0.re = in.e0.re * factor;
    tmp.e0.im = in.e0.im * factor;
    tmp.e1.re = in.e1.re * factor;
    tmp.e1.im = in.e1.im * factor;
    tmp.e2.re = in.e2.re * factor;
    tmp.e2.im = in.e2.im * factor;
    return tmp;
}

su3vec su3vec_times_complex(const su3vec in, const hmc_complex factor)
{
    su3vec tmp;
    tmp.e0.re = in.e0.re * factor.re - in.e0.im * factor.im;
    tmp.e0.im = in.e0.im * factor.re + in.e0.re * factor.im;
    tmp.e1.re = in.e1.re * factor.re - in.e1.im * factor.im;
    tmp.e1.im = in.e1.im * factor.re + in.e1.re * factor.im;
    tmp.e2.re = in.e2.re * factor.re - in.e2.im * factor.im;
    tmp.e2.im = in.e2.im * factor.re + in.e2.re * factor.im;
    return tmp;
}

su3vec su3vec_times_complex_conj(su3vec in, hmc_complex factor)
{
    su3vec tmp;
    tmp.e0.re = in.e0.re * factor.re + in.e0.im * factor.im;
    tmp.e0.im = in.e0.im * factor.re - in.e0.re * factor.im;
    tmp.e1.re = in.e1.re * factor.re + in.e1.im * factor.im;
    tmp.e1.im = in.e1.im * factor.re - in.e1.re * factor.im;
    tmp.e2.re = in.e2.re * factor.re + in.e2.im * factor.im;
    tmp.e2.im = in.e2.im * factor.re - in.e2.re * factor.im;
    return tmp;
}

su3vec su3matrix_times_su3vec(Matrixsu3 u, su3vec in)
{
    su3vec tmp;

    tmp.e0.re = u.e00.re * in.e0.re + u.e01.re * in.e1.re + u.e02.re * in.e2.re - u.e00.im * in.e0.im -
                u.e01.im * in.e1.im - u.e02.im * in.e2.im;
    tmp.e0.im = u.e00.re * in.e0.im + u.e01.re * in.e1.im + u.e02.re * in.e2.im + u.e00.im * in.e0.re +
                u.e01.im * in.e1.re + u.e02.im * in.e2.re;

    tmp.e1.re = u.e10.re * in.e0.re + u.e11.re * in.e1.re + u.e12.re * in.e2.re - u.e10.im * in.e0.im -
                u.e11.im * in.e1.im - u.e12.im * in.e2.im;
    tmp.e1.im = u.e10.re * in.e0.im + u.e11.re * in.e1.im + u.e12.re * in.e2.im + u.e10.im * in.e0.re +
                u.e11.im * in.e1.re + u.e12.im * in.e2.re;

    tmp.e2.re = u.e20.re * in.e0.re + u.e21.re * in.e1.re + u.e22.re * in.e2.re - u.e20.im * in.e0.im -
                u.e21.im * in.e1.im - u.e22.im * in.e2.im;
    tmp.e2.im = u.e20.re * in.e0.im + u.e21.re * in.e1.im + u.e22.re * in.e2.im + u.e20.im * in.e0.re +
                u.e21.im * in.e1.re + u.e22.im * in.e2.re;

    return tmp;
}

su3vec su3matrix_dagger_times_su3vec(Matrixsu3 u, su3vec in)
{
    su3vec tmp;

    tmp.e0.re = u.e00.re * in.e0.re + u.e10.re * in.e1.re + u.e20.re * in.e2.re + u.e00.im * in.e0.im +
                u.e10.im * in.e1.im + u.e20.im * in.e2.im;
    tmp.e0.im = u.e00.re * in.e0.im + u.e10.re * in.e1.im + u.e20.re * in.e2.im - u.e00.im * in.e0.re -
                u.e10.im * in.e1.re - u.e20.im * in.e2.re;

    tmp.e1.re = u.e01.re * in.e0.re + u.e11.re * in.e1.re + u.e21.re * in.e2.re + u.e01.im * in.e0.im +
                u.e11.im * in.e1.im + u.e21.im * in.e2.im;
    tmp.e1.im = u.e01.re * in.e0.im + u.e11.re * in.e1.im + u.e21.re * in.e2.im - u.e01.im * in.e0.re -
                u.e11.im * in.e1.re - u.e21.im * in.e2.re;

    tmp.e2.re = u.e02.re * in.e0.re + u.e12.re * in.e1.re + u.e22.re * in.e2.re + u.e02.im * in.e0.im +
                u.e12.im * in.e1.im + u.e22.im * in.e2.im;
    tmp.e2.im = u.e02.re * in.e0.im + u.e12.re * in.e1.im + u.e22.re * in.e2.im - u.e02.im * in.e0.re -
                u.e12.im * in.e1.re - u.e22.im * in.e2.re;

    return tmp;
}

su3vec su3vec_times_minusone(su3vec in)
{
    su3vec tmp;
    tmp.e0.re = -(in).e0.re;
    tmp.e0.im = -(in).e0.im;
    tmp.e1.re = -(in).e1.re;
    tmp.e1.im = -(in).e1.im;
    tmp.e2.re = -(in).e2.re;
    tmp.e2.im = -(in).e2.im;
    return tmp;
}

su3vec su3vec_acc(su3vec in1, su3vec in2)
{
    su3vec tmp;
    tmp.e0.re = in1.e0.re + in2.e0.re;
    tmp.e0.im = in1.e0.im + in2.e0.im;
    tmp.e1.re = in1.e1.re + in2.e1.re;
    tmp.e1.im = in1.e1.im + in2.e1.im;
    tmp.e2.re = in1.e2.re + in2.e2.re;
    tmp.e2.im = in1.e2.im + in2.e2.im;
    return tmp;
}

su3vec su3vec_acc_i(su3vec in1, su3vec in2)
{
    su3vec tmp;
    tmp.e0.re = in1.e0.re - in2.e0.im;
    tmp.e0.im = in1.e0.im + in2.e0.re;
    tmp.e1.re = in1.e1.re - in2.e1.im;
    tmp.e1.im = in1.e1.im + in2.e1.re;
    tmp.e2.re = in1.e2.re - in2.e2.im;
    tmp.e2.im = in1.e2.im + in2.e2.re;
    return tmp;
}

su3vec su3vec_dim(su3vec in1, su3vec in2)
{
    su3vec tmp;
    tmp.e0.re = in1.e0.re - in2.e0.re;
    tmp.e0.im = in1.e0.im - in2.e0.im;
    tmp.e1.re = in1.e1.re - in2.e1.re;
    tmp.e1.im = in1.e1.im - in2.e1.im;
    tmp.e2.re = in1.e2.re - in2.e2.re;
    tmp.e2.im = in1.e2.im - in2.e2.im;
    return tmp;
}

su3vec su3vec_dim_i(su3vec in1, su3vec in2)
{
    su3vec tmp;
    tmp.e0.re = in1.e0.re + in2.e0.im;
    tmp.e0.im = in1.e0.im - in2.e0.re;
    tmp.e1.re = in1.e1.re + in2.e1.im;
    tmp.e1.im = in1.e1.im - in2.e1.re;
    tmp.e2.re = in1.e2.re + in2.e2.im;
    tmp.e2.im = in1.e2.im - in2.e2.re;
    return tmp;
}

su3vec su3vec_acc_acc(su3vec in1, su3vec in2, su3vec in3)
{
    su3vec tmp;
    tmp.e0.re = in1.e0.re + in2.e0.re + in3.e0.re;
    tmp.e0.im = in1.e0.im + in2.e0.im + in3.e0.im;
    tmp.e1.re = in1.e1.re + in2.e1.re + in3.e1.re;
    tmp.e1.im = in1.e1.im + in2.e1.im + in3.e1.im;
    tmp.e2.re = in1.e2.re + in2.e2.re + in3.e2.re;
    tmp.e2.im = in1.e2.im + in2.e2.im + in3.e2.im;
    return tmp;
}

// calculates the Dirac-Trace of the matrix resulting from multiplying U*V^dagger =  u*v^dagger + w*x^dagger, where u,
// v, w, x are SU(3)-vectors (using spinprojection). The result is a 3x3-matrix
Matrix3x3 tr_v_times_u_dagger(su3vec u, su3vec v, su3vec w, su3vec x)
{
    Matrix3x3 tmp;
    tmp.e00.re = (u).e0.re * (v).e0.re + (u).e0.im * (v).e0.im + (w).e0.re * (x).e0.re + (w).e0.im * (x).e0.im;
    tmp.e00.im = (u).e0.im * (v).e0.re - (u).e0.re * (v).e0.im + (w).e0.im * (x).e0.re - (w).e0.re * (x).e0.im;
    tmp.e01.re = (u).e0.re * (v).e1.re + (u).e0.im * (v).e1.im + (w).e0.re * (x).e1.re + (w).e0.im * (x).e1.im;
    tmp.e01.im = (u).e0.im * (v).e1.re - (u).e0.re * (v).e1.im + (w).e0.im * (x).e1.re - (w).e0.re * (x).e1.im;
    tmp.e02.re = (u).e0.re * (v).e2.re + (u).e0.im * (v).e2.im + (w).e0.re * (x).e2.re + (w).e0.im * (x).e2.im;
    tmp.e02.im = (u).e0.im * (v).e2.re - (u).e0.re * (v).e2.im + (w).e0.im * (x).e2.re - (w).e0.re * (x).e2.im;
    tmp.e10.re = (u).e1.re * (v).e0.re + (u).e1.im * (v).e0.im + (w).e1.re * (x).e0.re + (w).e1.im * (x).e0.im;
    tmp.e10.im = (u).e1.im * (v).e0.re - (u).e1.re * (v).e0.im + (w).e1.im * (x).e0.re - (w).e1.re * (x).e0.im;
    tmp.e11.re = (u).e1.re * (v).e1.re + (u).e1.im * (v).e1.im + (w).e1.re * (x).e1.re + (w).e1.im * (x).e1.im;
    tmp.e11.im = (u).e1.im * (v).e1.re - (u).e1.re * (v).e1.im + (w).e1.im * (x).e1.re - (w).e1.re * (x).e1.im;
    tmp.e12.re = (u).e1.re * (v).e2.re + (u).e1.im * (v).e2.im + (w).e1.re * (x).e2.re + (w).e1.im * (x).e2.im;
    tmp.e12.im = (u).e1.im * (v).e2.re - (u).e1.re * (v).e2.im + (w).e1.im * (x).e2.re - (w).e1.re * (x).e2.im;
    tmp.e20.re = (u).e2.re * (v).e0.re + (u).e2.im * (v).e0.im + (w).e2.re * (x).e0.re + (w).e2.im * (x).e0.im;
    tmp.e20.im = (u).e2.im * (v).e0.re - (u).e2.re * (v).e0.im + (w).e2.im * (x).e0.re - (w).e2.re * (x).e0.im;
    tmp.e21.re = (u).e2.re * (v).e1.re + (u).e2.im * (v).e1.im + (w).e2.re * (x).e1.re + (w).e2.im * (x).e1.im;
    tmp.e21.im = (u).e2.im * (v).e1.re - (u).e2.re * (v).e1.im + (w).e2.im * (x).e1.re - (w).e2.re * (x).e1.im;
    tmp.e22.re = (u).e2.re * (v).e2.re + (u).e2.im * (v).e2.im + (w).e2.re * (x).e2.re + (w).e2.im * (x).e2.im;
    tmp.e22.im = (u).e2.im * (v).e2.re - (u).e2.re * (v).e2.im + (w).e2.im * (x).e2.re - (w).e2.re * (x).e2.im;
    return tmp;
}

// This function makes the direct product of the su3 vectors u and v, returning the Matrix3x3 u*v^dagger
Matrix3x3 u_times_v_dagger(su3vec u, su3vec v)
{
    Matrix3x3 out;

    out.e00.re = u.e0.re * v.e0.re + u.e0.im * v.e0.im;
    out.e00.im = u.e0.im * v.e0.re - u.e0.re * v.e0.im;
    out.e01.re = u.e0.re * v.e1.re + u.e0.im * v.e1.im;
    out.e01.im = u.e0.im * v.e1.re - u.e0.re * v.e1.im;
    out.e02.re = u.e0.re * v.e2.re + u.e0.im * v.e2.im;
    out.e02.im = u.e0.im * v.e2.re - u.e0.re * v.e2.im;
    out.e10.re = u.e1.re * v.e0.re + u.e1.im * v.e0.im;
    out.e10.im = u.e1.im * v.e0.re - u.e1.re * v.e0.im;
    out.e11.re = u.e1.re * v.e1.re + u.e1.im * v.e1.im;
    out.e11.im = u.e1.im * v.e1.re - u.e1.re * v.e1.im;
    out.e12.re = u.e1.re * v.e2.re + u.e1.im * v.e2.im;
    out.e12.im = u.e1.im * v.e2.re - u.e1.re * v.e2.im;
    out.e20.re = u.e2.re * v.e0.re + u.e2.im * v.e0.im;
    out.e20.im = u.e2.im * v.e0.re - u.e2.re * v.e0.im;
    out.e21.re = u.e2.re * v.e1.re + u.e2.im * v.e1.im;
    out.e21.im = u.e2.im * v.e1.re - u.e2.re * v.e1.im;
    out.e22.re = u.e2.re * v.e2.re + u.e2.im * v.e2.im;
    out.e22.im = u.e2.im * v.e2.re - u.e2.re * v.e2.im;

    return out;
}
