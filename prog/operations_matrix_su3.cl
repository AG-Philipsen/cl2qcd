/** @file
 * Device code implementing SU(3)matrices with and without reconstruct 12
 */
//operations_matrix_su3.cl

#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(const Matrixsu3 p, const int ncomp)
{
	hmc_complex out;
	hmc_complex tmp;

	switch (ncomp) {

		case 0:
			tmp = complexmult (p.e01, p.e12);
			out = complexmult (p.e02, p.e11);
			break;

		case 1:
			tmp = complexmult (p.e02, p.e10);
			out = complexmult (p.e00, p.e12);
			break;
		case 2:
			tmp = complexmult (p.e00, p.e11);
			out = complexmult (p.e01, p.e10);
			break;
	}

	out = complexsubtract (tmp, out);
	return complexconj (out);
}
#endif

void print_matrixsu3(Matrixsu3 in)
{
#ifdef _RECONSTRUCT_TWELVE_
	hmc_float in20_re = reconstruct_su3(in, 0).re;
	hmc_float in20_im = reconstruct_su3(in, 0).im;
	hmc_float in21_re = reconstruct_su3(in, 1).re;
	hmc_float in21_im = reconstruct_su3(in, 1).im;
	hmc_float in22_re = reconstruct_su3(in, 2).re;
	hmc_float in22_im = reconstruct_su3(in, 2).im;

	printf("(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n",
	       in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im,
	       in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im,
	       in20_re, in20_im, in21_re, in21_im, in22_re, in22_im);
#else
	printf("(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n",
	       in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im,
	       in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im,
	       in.e20.re, in.e20.im, in.e21.re, in.e21.im, in.e22.re, in.e22.im);
#endif
	printf("\n");
}



Matrixsu3 get_matrixsu3( __global ocl_s_gaugefield * field, const int spacepos, const int timepos, const int mu)
{
	Matrixsu3  out;
	out = field [get_global_link_pos(mu, spacepos, timepos)];


	return out;
}

void put_matrixsu3(__global ocl_s_gaugefield * field, const Matrixsu3 in, const int spacepos, const int timepos, const int mu)
{
	field [get_global_link_pos(mu, spacepos, timepos)] = in;
}

Matrixsu3 copy_matrixsu3(const Matrixsu3 in)
{
	Matrixsu3 out;

	out.e00.re = in.e00.re;
	out.e00.im = in.e00.im;
	out.e01.re = in.e01.re;
	out.e01.im = in.e01.im;
	out.e02.re = in.e02.re;
	out.e02.im = in.e02.im;

	out.e10.re = in.e10.re;
	out.e10.im = in.e10.im;
	out.e11.re = in.e11.re;
	out.e11.im = in.e11.im;
	out.e12.re = in.e12.re;
	out.e12.im = in.e12.im;

#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re = in.e20.re;
	out.e20.im = in.e20.im;
	out.e21.re = in.e21.re;
	out.e21.im = in.e21.im;
	out.e22.re = in.e22.re;
	out.e22.im = in.e22.im;
#endif
	return out;
}

Matrixsu3 unit_matrixsu3()
{
	Matrixsu3 out;
	out.e00.re = 1.;
	out.e00.im = 0.;
	out.e01.re = 0.;
	out.e01.im = 0.;
	out.e02.re = 0.;
	out.e02.im = 0.;

	out.e10.re = 0.;
	out.e10.im = 0.;
	out.e11.re = 1.;
	out.e11.im = 0.;
	out.e12.re = 0.;
	out.e12.im = 0.;

#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 1.;
	out.e22.im = 0.;
#endif
	return out;
}


Matrixsu3 zero_matrixsu3()
{
	Matrixsu3 out;
	out.e00.re = 0.;
	out.e00.im = 0.;
	out.e01.re = 0.;
	out.e01.im = 0.;
	out.e02.re = 0.;
	out.e02.im = 0.;

	out.e10.re = 0.;
	out.e10.im = 0.;
	out.e11.re = 0.;
	out.e11.im = 0.;
	out.e12.re = 0.;
	out.e12.im = 0.;

#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 0.;
	out.e22.im = 0.;
#endif

	return out;
}

Matrixsu3 multiply_matrixsu3(const Matrixsu3 p, const Matrixsu3 q)
{
	Matrixsu3 out;
#ifdef _RECONSTRUCT_TWELVE_
	hmc_float p20_re = reconstruct_su3(p, 0).re;
	hmc_float p20_im = reconstruct_su3(p, 0).im;
	hmc_float p21_re = reconstruct_su3(p, 1).re;
	hmc_float p21_im = reconstruct_su3(p, 1).im;
	hmc_float p22_re = reconstruct_su3(p, 2).re;
	hmc_float p22_im = reconstruct_su3(p, 2).im;

	hmc_float q20_re = reconstruct_su3(q, 0).re;
	hmc_float q20_im = reconstruct_su3(q, 0).im;
	hmc_float q21_re = reconstruct_su3(q, 1).re;
	hmc_float q21_im = reconstruct_su3(q, 1).im;
	hmc_float q22_re = reconstruct_su3(q, 2).re;
	hmc_float q22_im = reconstruct_su3(q, 2).im;

	out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e10.re + p.e02.re * q20_re
	             - p.e00.im * q.e00.im - p.e01.im * q.e10.im - p.e02.im * q20_im;
	out.e00.im = p.e00.re * q.e00.im + p.e01.re * q.e10.im + p.e02.re * q20_im
	             + p.e00.im * q.e00.re + p.e01.im * q.e10.re + p.e02.im * q20_re;

	out.e01.re = p.e00.re * q.e01.re + p.e01.re * q.e11.re + p.e02.re * q21_re
	             - p.e00.im * q.e01.im - p.e01.im * q.e11.im - p.e02.im * q21_im;
	out.e01.im = p.e00.re * q.e01.im + p.e01.re * q.e11.im + p.e02.re * q21_im
	             + p.e00.im * q.e01.re + p.e01.im * q.e11.re + p.e02.im * q21_re;

	out.e02.re = p.e00.re * q.e02.re + p.e01.re * q.e12.re + p.e02.re * q22_re
	             - p.e00.im * q.e02.im - p.e01.im * q.e12.im - p.e02.im * q22_im;
	out.e02.im = p.e00.re * q.e02.im + p.e01.re * q.e12.im + p.e02.re * q22_im
	             + p.e00.im * q.e02.re + p.e01.im * q.e12.re + p.e02.im * q22_re;

	out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e10.re + p.e12.re * q20_re
	             - p.e10.im * q.e00.im - p.e11.im * q.e10.im - p.e12.im * q20_im;
	out.e10.im = p.e10.re * q.e00.im + p.e11.re * q.e10.im + p.e12.re * q20_im
	             + p.e10.im * q.e00.re + p.e11.im * q.e10.re + p.e12.im * q20_re;

	out.e11.re = p.e10.re * q.e01.re + p.e11.re * q.e11.re + p.e12.re * q21_re
	             - p.e10.im * q.e01.im - p.e11.im * q.e11.im - p.e12.im * q21_im;
	out.e11.im = p.e10.re * q.e01.im + p.e11.re * q.e11.im + p.e12.re * q21_im
	             + p.e10.im * q.e01.re + p.e11.im * q.e11.re + p.e12.im * q21_re;

	out.e12.re = p.e10.re * q.e02.re + p.e11.re * q.e12.re + p.e12.re * q22_re
	             - p.e10.im * q.e02.im - p.e11.im * q.e12.im - p.e12.im * q22_im;
	out.e12.im = p.e10.re * q.e02.im + p.e11.re * q.e12.im + p.e12.re * q22_im
	             + p.e10.im * q.e02.re + p.e11.im * q.e12.re + p.e12.im * q22_re;
#else

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

#endif

	return out;
}

Matrixsu3 multiply_matrixsu3_dagger(const Matrixsu3 p, const Matrixsu3 q)
{
	Matrixsu3 out;
#ifdef _RECONSTRUCT_TWELVE_
	hmc_float p20_re = reconstruct_su3(p, 0).re;
	hmc_float p20_im = reconstruct_su3(p, 0).im;
	hmc_float p21_re = reconstruct_su3(p, 1).re;
	hmc_float p21_im = reconstruct_su3(p, 1).im;
	hmc_float p22_re = reconstruct_su3(p, 2).re;
	hmc_float p22_im = reconstruct_su3(p, 2).im;

	hmc_float q20_re = reconstruct_su3(q, 0).re;
	hmc_float q20_im = reconstruct_su3(q, 0).im;
	hmc_float q21_re = reconstruct_su3(q, 1).re;
	hmc_float q21_im = reconstruct_su3(q, 1).im;
	hmc_float q22_re = reconstruct_su3(q, 2).re;
	hmc_float q22_im = reconstruct_su3(q, 2).im;

	out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
	             + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
	out.e00.im = -p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
	             + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

	out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
	             + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
	out.e01.im = -p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
	             + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;

	out.e02.re = p.e00.re * q20_re   + p.e01.re * q21_re   + p.e02.re * q22_re
	             + p.e00.im * q20_im   + p.e01.im * q21_im   + p.e02.im * q22_im;
	out.e02.im = -p.e00.re * q20_im   - p.e01.re * q21_im   - p.e02.re * q22_im
	             + p.e00.im * q20_re   + p.e01.im * q21_re   + p.e02.im * q22_re;

	out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
	             + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
	out.e10.im = -p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
	             + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;

	out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
	             + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
	out.e11.im = -p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
	             + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;

	out.e12.re = p.e10.re * q20_re   + p.e11.re * q21_re   + p.e12.re * q22_re
	             + p.e10.im * q20_im   + p.e11.im * q21_im   + p.e12.im * q22_im;
	out.e12.im = -p.e10.re * q20_im   - p.e11.re * q21_im   - p.e12.re * q22_im
	             + p.e10.im * q20_re   + p.e11.im * q21_re   + p.e12.im * q22_re;
#else

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

#endif
	return out;
}

Matrixsu3 adjoint_matrixsu3(const Matrixsu3 p)
{
	Matrixsu3 out;
	out.e00.re = p.e00.re;
	out.e00.im = - p.e00.im;
	out.e01.re = p.e10.re;
	out.e01.im = - p.e10.im;

	out.e10.re = p.e01.re;
	out.e10.im = - p.e01.im;
	out.e11.re = p.e11.re;
	out.e11.im = - p.e11.im;


#ifdef _RECONSTRUCT_TWELVE_
	out.e02.re = reconstruct_su3(p, 0).re;
	out.e02.im = - reconstruct_su3(p, 0).im;

	out.e12.re = reconstruct_su3(p, 1).re;
	out.e12.im = - reconstruct_su3(p, 1).im;
#else
	out.e02.re = p.e20.re;
	out.e02.im = - p.e20.im;

	out.e12.re = p.e21.re;
	out.e12.im = - p.e21.im;

	out.e20.re = p.e02.re;
	out.e20.im = - p.e02.im;
	out.e21.re = p.e12.re;
	out.e21.im = - p.e12.im;
	out.e22.re = p.e22.re;
	out.e22.im = - p.e22.im;
#endif
	return out;
}

hmc_complex trace_matrixsu3(const Matrixsu3 p)
{
	hmc_complex out;
	out.re = p.e00.re;
	out.im = p.e00.im;
	out.re += p.e11.re;
	out.im += p.e11.im;
#ifdef _RECONSTRUCT_TWELVE_
	out.re += reconstruct_su3(p, 2).re;
	out.im += reconstruct_su3(p, 2).im;
#else
	out.re += p.e22.re;
	out.im += p.e22.im;
#endif
	return out;
}

hmc_complex det_matrixsu3(const Matrixsu3 p)
{

	hmc_complex out;

#ifdef _RECONSTRUCT_TWELVE_
	hmc_complex c0 = reconstruct_su3(p, 0);
	hmc_complex c1 = reconstruct_su3(p, 1);
	hmc_complex c2 = reconstruct_su3(p, 2);


	hmc_complex det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	tmp1 = complexmult( p.e11, c2);
	det1 = complexmult( p.e00, tmp1);
	tmp2 = complexmult( p.e12, c0);
	det2 = complexmult( p.e01, tmp2);
	tmp3 = complexmult( p.e10, c1);
	det3 = complexmult( p.e02, tmp3);
	tmp4 = complexmult( p.e11, c0 );
	det4 = complexmult( p.e02, tmp4);
	tmp5 = complexmult( p.e10, c2 );
	det5 = complexmult( p.e01, tmp5);
	tmp6 = complexmult( p.e12, c1 );
	det6 = complexmult( p.e00, tmp6);
#else
	hmc_complex det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	tmp1 = complexmult( p.e11, p.e22);
	det1 = complexmult( p.e00, tmp1);
	tmp2 = complexmult( p.e12, p.e20);
	det2 = complexmult( p.e01, tmp2);
	tmp3 = complexmult( p.e10, p.e21);
	det3 = complexmult( p.e02, tmp3);
	tmp4 = complexmult( p.e11, p.e20 );
	det4 = complexmult( p.e02, tmp4);
	tmp5 = complexmult( p.e10, p.e22 );
	det5 = complexmult( p.e01, tmp5);
	tmp6 = complexmult( p.e12, p.e21 );
	det6 = complexmult( p.e00, tmp6);
#endif
	out.re = det1.re + det2.re + det3.re - det4.re - det5.re - det6.re;
	out.im = det1.im + det2.im + det3.im - det4.im - det5.im - det6.im;

	return out;
}

//CP: tested version that recreates the matrices made by tmlqcd
/** @todo recheck the factor 0.5 (or F_1_2) that has been deleted here */

//this build a su3-matrix from an algebraelement!!
Matrixsu3 build_su3_from_ae(ae in)
{
	Matrixsu3 v;

	v.e00.re = 0.0;
	v.e00.im = (in.e7 * F_1_S3 + in.e2);
	v.e01.re = in.e1;
	v.e01.im = in.e0;
	v.e02.re = in.e4;
	v.e02.im = in.e3;
	v.e10.re = -in.e1;;
	v.e10.im = in.e0;;
	v.e11.re = 0.0;
	v.e11.im = (in.e7 * F_1_S3 - in.e2);
	v.e12.re = in.e6;
	v.e12.im = in.e5;
#ifndef _RECONSTRUCT_TWELVE_
	v.e20.re = -in.e4;
	v.e20.im = in.e3;
	v.e21.re = -in.e6;
	v.e21.im = in.e5;
	v.e22.re = 0.0;
	v.e22.im = -2.*in.e7 * F_1_S3;
#endif
	return v;
}

//this build a su3-matrix from an algebraelement and in addition multiplies it by a real number!!
Matrixsu3 build_su3_from_ae_times_real(ae in, hmc_float eps)
{
	Matrixsu3 v;

	v.e00.re = 0.0;
	v.e00.im = eps * (in.e7 * F_1_S3 + in.e2);
	v.e01.re = eps * in.e1;
	v.e01.im = eps * in.e0;
	v.e02.re = eps * in.e4;
	v.e02.im = eps * in.e3;
	v.e10.re = -eps * in.e1;;
	v.e10.im = eps * in.e0;;
	v.e11.re = 0.0;
	v.e11.im = eps * (in.e7 * F_1_S3 - in.e2);
	v.e12.re = eps * in.e6;
	v.e12.im = eps * in.e5;
#ifndef _RECONSTRUCT_TWELVE_
	v.e20.re = -eps * in.e4;
	v.e20.im = eps * in.e3;
	v.e21.re = -eps * in.e6;
	v.e21.im = eps * in.e5;
	v.e22.re = 0.0;
	v.e22.im = -eps * 2.*in.e7 * F_1_S3;
#endif

	return v;
}

Matrixsu3 build_su3matrix_by_exponentiation(ae inn, hmc_float epsilon)
{
	//CP: this is the old version

//	Matrixsu3 tmp;
//	hmc_float in [8];
//	in[0] = inn.e0;
//	in[1] = inn.e1;
//	in[2] = inn.e2;
//	in[3] = inn.e3;
//	in[4] = inn.e4;
//	in[5] = inn.e5;
//	in[6] = inn.e6;
//	in[7] = inn.e7;
//
//	//this is the case where one actually evaluates -as matrices- many orders of exp(i*e*P)=1+i*e*P+(1/2)(i*e*P)^2 + ...
//	//1. create mat = (i epsilon)/2 p_i lambda_i    with lambda=gellmann matrix
//	Matrix3x3 eMat;
//	hmc_float halfeps = epsilon;//*F_1_2;
//	eMat.e00.re = 0.0;
//	eMat.e00.im = halfeps * (in[7] * F_1_S3 + in[2]);
//	eMat.e01.re = halfeps * in[1];
//	eMat.e01.im = halfeps * in[0];
//	eMat.e02.re = halfeps * in[4];
//	eMat.e02.im = halfeps * in[3];
//	eMat.e10.re = -halfeps * in[1];;
//	eMat.e10.im = halfeps * in[0];;
//	eMat.e11.re = 0.0;
//	eMat.e11.im = halfeps * (in[7] * F_1_S3 - in[2]);
//	eMat.e12.re = halfeps * in[6];
//	eMat.e12.im = halfeps * in[5];
//	eMat.e20.re = -halfeps * in[4];
//	eMat.e20.im = halfeps * in[3];
//	eMat.e21.re = -halfeps * in[6];
//	eMat.e21.im = halfeps * in[5];
//	eMat.e22.re = 0.0;
//	eMat.e22.im = -halfeps * 2.*in[7] * F_1_S3;
//
//	//2. start with the exp(...) expansion by using standard 3x3 sum and multiplication
//	Matrix3x3 eRes, eCurPower, eNextPower, eLastResult;
//	hmc_float eAccuracyCheck;
//	eRes = identity_matrix3x3();
//	eCurPower = identity_matrix3x3();
//	hmc_float eCoefficient = 1.0;
//	int factorial = 1;
//	for(int power = 1; power < 50; power++) {
//		eNextPower = multiply_matrix3x3(eMat, eCurPower);
//		eCurPower = eNextPower;
//		//CP: the calc of factorial can be simplified by just doing factorial*=power!!
//		for(int i = power; i > 1; i--) factorial *= i;
//		eCurPower = multiply_matrix3x3_by_real(eCurPower, (1. / (1.*factorial)));
//		factorial = 1;
//		eLastResult = eRes;
//		eRes = add_matrix3x3(eRes, eCurPower);
//		eAccuracyCheck = absoluteDifference_matrix3x3(eRes, eLastResult);
//		if(eAccuracyCheck < 1e-15) {
//			break;
//		}
//	}
//	//3. here I have the exponentiated matrix in 3x3 generic form (eRes), project it
//#ifdef _RECONSTRUCT_TWELVE_
//	tmp.e0 = eRes.e00;
//	tmp.e1 = eRes.e01;
//	tmp.e2 = eRes.e02;
//	tmp.e3 = eRes.e10;
//	tmp.e4 = eRes.e11;
//	tmp.e5 = eRes.e12;
//#else
//	tmp.e00 = eRes.e00;
//	tmp.e01 = eRes.e01;
//	tmp.e02 = eRes.e02;
//	tmp.e10 = eRes.e10;
//	tmp.e11 = eRes.e11;
//	tmp.e12 = eRes.e12;
//	tmp.e20 = eRes.e20;
//	tmp.e21 = eRes.e21;
//	tmp.e22 = eRes.e22;
//#endif // _RECONSTRUCT_TWELVE_
//	project_su3(tmp);
//
//	return tmp;
//

	//this is the method taken from tmqlcd. It is a cut-off series.
	int i;
	Matrixsu3 v, v2, vr;
	hmc_float fac, r;
	hmc_float a, b;
	hmc_complex a0, a1, a2, a1p;

	//original in tmlqcd: _make_su3(v,p);
	//Here, the stepsize factor is multiplied in right away!!
	//Also, a factor of 0.5 is taken out to fit the different trace-definition from tmqlcd
	hmc_float halfeps = epsilon;// *F_1_2;
	v = build_su3_from_ae_times_real(inn, halfeps);

	// calculates v^2
	v2 = multiply_matrixsu3(v, v);
#ifdef _RECONSTRUCT_TWELVE_
	hmc_complex v_e20 = reconstruct_su3(v, 0);
	hmc_complex v_e21 = reconstruct_su3(v, 1);
	hmc_complex v_e22 = reconstruct_su3(v, 2);

	hmc_complex v2_e20 = reconstruct_su3(v2, 0);
	hmc_complex v2_e21 = reconstruct_su3(v2, 1);
	hmc_complex v2_e22 = reconstruct_su3(v2, 2);

	a = 0.5 * (v2.e00.re + v2.e11.re + v2_e22.re);

	b = 0.33333333333333333 *
	    (v.e00.re * v2.e00.im + v.e00.im * v2.e00.re
	     + v.e01.re * v2.e10.im + v.e01.im * v2.e10.re
	     + v.e02.re * v2_e20.im + v.e02.im * v2_e20.re
	     + v.e10.re * v2.e01.im + v.e10.im * v2.e01.re
	     + v.e11.re * v2.e11.im + v.e11.im * v2.e11.re
	     + v.e12.re * v2_e21.im + v.e12.im * v2_e21.re
	     + v_e20.re * v2.e02.im + v_e20.im * v2.e02.re
	     + v_e21.re * v2.e12.im + v_e21.im * v2.e12.re
	     + v_e22.re * v2_e22.im + v_e22.im * v2_e22.re  );

#else
	a = 0.5 * (v2.e00.re + v2.e11.re + v2.e22.re);
	// 1/3 imaginary part of tr v*v2
	b = 0.33333333333333333 *
	    (v.e00.re * v2.e00.im + v.e00.im * v2.e00.re
	     + v.e01.re * v2.e10.im + v.e01.im * v2.e10.re
	     + v.e02.re * v2.e20.im + v.e02.im * v2.e20.re
	     + v.e10.re * v2.e01.im + v.e10.im * v2.e01.re
	     + v.e11.re * v2.e11.im + v.e11.im * v2.e11.re
	     + v.e12.re * v2.e21.im + v.e12.im * v2.e21.re
	     + v.e20.re * v2.e02.im + v.e20.im * v2.e02.re
	     + v.e21.re * v2.e12.im + v.e21.im * v2.e12.re
	     + v.e22.re * v2.e22.im + v.e22.im * v2.e22.re  );
#endif
	a0.re = 0.16059043836821615e-9;  //  1/13!
	a0.im = 0.0;
	a1.re = 0.11470745597729725e-10; //  1/14!
	a1.im = 0.0;
	a2.re = 0.76471637318198165e-12; //  1/15!
	a2.im = 0.0;
	fac = 0.20876756987868099e-8;    //  1/12!
	r = 12.0;
	for(i = 3; i <= 15; i++) {
		a1p.re = a0.re + a * a2.re;
		a1p.im = a0.im + a * a2.im;
		a0.re = fac - b * a2.im;
		a0.im =     + b * a2.re;
		a2.re = a1.re;
		a2.im = a1.im;
		a1.re = a1p.re;
		a1.im = a1p.im;
		fac *= r;
		r -= 1.0;
	}
	// vr = a0 + a1*v + a2*v2
	vr.e00.re = a0.re + a1.re * v.e00.re - a1.im * v.e00.im + a2.re * v2.e00.re - a2.im * v2.e00.im;
	vr.e00.im = a0.im + a1.re * v.e00.im + a1.im * v.e00.re + a2.re * v2.e00.im + a2.im * v2.e00.re;
	vr.e01.re =         a1.re * v.e01.re - a1.im * v.e01.im + a2.re * v2.e01.re - a2.im * v2.e01.im;
	vr.e01.im =         a1.re * v.e01.im + a1.im * v.e01.re + a2.re * v2.e01.im + a2.im * v2.e01.re;
	vr.e02.re =         a1.re * v.e02.re - a1.im * v.e02.im + a2.re * v2.e02.re - a2.im * v2.e02.im;
	vr.e02.im =         a1.re * v.e02.im + a1.im * v.e02.re + a2.re * v2.e02.im + a2.im * v2.e02.re;
	vr.e10.re =         a1.re * v.e10.re - a1.im * v.e10.im + a2.re * v2.e10.re - a2.im * v2.e10.im;
	vr.e10.im =         a1.re * v.e10.im + a1.im * v.e10.re + a2.re * v2.e10.im + a2.im * v2.e10.re;
	vr.e11.re = a0.re + a1.re * v.e11.re - a1.im * v.e11.im + a2.re * v2.e11.re - a2.im * v2.e11.im;
	vr.e11.im = a0.im + a1.re * v.e11.im + a1.im * v.e11.re + a2.re * v2.e11.im + a2.im * v2.e11.re;
	vr.e12.re =         a1.re * v.e12.re - a1.im * v.e12.im + a2.re * v2.e12.re - a2.im * v2.e12.im;
	vr.e12.im =         a1.re * v.e12.im + a1.im * v.e12.re + a2.re * v2.e12.im + a2.im * v2.e12.re;
#ifndef _RECONSTRUCT_TWELVE_
	vr.e20.re =         a1.re * v.e20.re - a1.im * v.e20.im + a2.re * v2.e20.re - a2.im * v2.e20.im;
	vr.e20.im =         a1.re * v.e20.im + a1.im * v.e20.re + a2.re * v2.e20.im + a2.im * v2.e20.re;
	vr.e21.re =         a1.re * v.e21.re - a1.im * v.e21.im + a2.re * v2.e21.re - a2.im * v2.e21.im;
	vr.e21.im =         a1.re * v.e21.im + a1.im * v.e21.re + a2.re * v2.e21.im + a2.im * v2.e21.re;
	vr.e22.re = a0.re + a1.re * v.e22.re - a1.im * v.e22.im + a2.re * v2.e22.re - a2.im * v2.e22.im;
	vr.e22.im = a0.im + a1.re * v.e22.im + a1.im * v.e22.re + a2.re * v2.e22.im + a2.im * v2.e22.re;
#endif
	return vr;

}

//scale a su3 matrix by a real factor
Matrixsu3 multiply_matrixsu3_by_real (Matrixsu3 in, hmc_float factor)
{
	Matrixsu3 out = in;
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
#ifndef _RECONSTRUCT_TWELVE_
	out.e20.re *= factor;
	out.e20.im *= factor;
	out.e21.re *= factor;
	out.e21.im *= factor;
	out.e22.re *= factor;
	out.e22.im *= factor;
#endif

	return out;
}

Matrixsu3 multiply_matrixsu3_by_complex (Matrixsu3 in, hmc_complex factor)
{
	Matrixsu3 out;
	out.e00 = complexmult(in.e00, factor);
	out.e01 = complexmult(in.e01, factor);
	out.e02 = complexmult(in.e02, factor);
	out.e10 = complexmult(in.e10, factor);
	out.e11 = complexmult(in.e11, factor);
	out.e12 = complexmult(in.e12, factor);
#ifndef _RECONSTRUCT_TWELVE_
	out.e20 = complexmult(in.e20, factor);
	out.e21 = complexmult(in.e21, factor);
	out.e22 = complexmult(in.e22, factor);
#endif
	return out;
}

