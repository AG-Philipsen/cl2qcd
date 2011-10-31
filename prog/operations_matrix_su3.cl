/** @file
 * Device code implementing SU(3)matrices with and without reconstruct 12
 */
//operations_matrix_su3.cl

#ifdef cl_amd_printf
void print_matrixsu3(Matrixsu3 in)
{
	printf("(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n",
	       in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im,
	       in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im,
	       in.e20.re, in.e20.im, in.e21.re, in.e21.im, in.e22.re, in.e22.im);
	printf("\n");
}
#endif



Matrixsu3 get_matrixsu3( __global ocl_s_gaugefield * field, const int spacepos, const int timepos, const int mu)
{
	return field [get_global_link_pos(mu, spacepos, timepos)];
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

	out.e20.re = in.e20.re;
	out.e20.im = in.e20.im;
	out.e21.re = in.e21.re;
	out.e21.im = in.e21.im;
	out.e22.re = in.e22.re;
	out.e22.im = in.e22.im;

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

	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 1.;
	out.e22.im = 0.;

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

	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 0.;
	out.e22.im = 0.;

	return out;
}

Matrixsu3 multiply_matrixsu3(const Matrixsu3 p, const Matrixsu3 q)
{
	Matrixsu3 out;

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

Matrixsu3 multiply_matrixsu3_dagger(const Matrixsu3 p, const Matrixsu3 q)
{
	Matrixsu3 out;

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

	return out;
}

hmc_complex trace_matrixsu3(const Matrixsu3 p)
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

hmc_complex det_matrixsu3(const Matrixsu3 p)
{

	hmc_complex out;

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
	v.e20.re = -in.e4;
	v.e20.im = in.e3;
	v.e21.re = -in.e6;
	v.e21.im = in.e5;
	v.e22.re = 0.0;
	v.e22.im = -2.*in.e7 * F_1_S3;
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
	v.e20.re = -eps * in.e4;
	v.e20.im = eps * in.e3;
	v.e21.re = -eps * in.e6;
	v.e21.im = eps * in.e5;
	v.e22.re = 0.0;
	v.e22.im = -eps * 2.*in.e7 * F_1_S3;

	return v;
}

Matrixsu3 build_su3matrix_by_exponentiation(ae inn, hmc_float epsilon)
{
	//this is the method taken from tmqlcd. It is a cut-off series.
	//original in tmlqcd: _make_su3(v,p);
	//Here, the stepsize factor is multiplied in right away!!
	//Also, a factor of 0.5 is taken out to fit the different trace-definition from tmqlcd
	hmc_float halfeps = epsilon;// *F_1_2;
	const Matrixsu3 v = build_su3_from_ae_times_real(inn, halfeps);

	// calculates v^2
	const Matrixsu3 v2 = multiply_matrixsu3(v, v);
	const hmc_float a = 0.5 * (v2.e00.re + v2.e11.re + v2.e22.re);
	// 1/3 imaginary part of tr v*v2
	const hmc_float b = 0.33333333333333333 *
	                    (v.e00.re * v2.e00.im + v.e00.im * v2.e00.re
	                     + v.e01.re * v2.e10.im + v.e01.im * v2.e10.re
	                     + v.e02.re * v2.e20.im + v.e02.im * v2.e20.re
	                     + v.e10.re * v2.e01.im + v.e10.im * v2.e01.re
	                     + v.e11.re * v2.e11.im + v.e11.im * v2.e11.re
	                     + v.e12.re * v2.e21.im + v.e12.im * v2.e21.re
	                     + v.e20.re * v2.e02.im + v.e20.im * v2.e02.re
	                     + v.e21.re * v2.e12.im + v.e21.im * v2.e12.re
	                     + v.e22.re * v2.e22.im + v.e22.im * v2.e22.re  );
	hmc_complex a0, a1, a2;
	a0.re = 0.16059043836821615e-9;  //  1/13!
	a0.im = 0.0;
	a1.re = 0.11470745597729725e-10; //  1/14!
	a1.im = 0.0;
	a2.re = 0.76471637318198165e-12; //  1/15!
	a2.im = 0.0;
	hmc_float fac = 0.20876756987868099e-8;
	uint r = 12;
	/*
	 * This pragma unroll is not a performance optimization but a required workaround for a bug in APP 2.5. It seems
	 * without it the GPU performs a wrong number of loop iterations.
	 * See http://code.compeng.uni-frankfurt.de/issues/211 and http://forums.amd.com/forum/messageview.cfm?catid=390&threadid=155815 .
	 */
#pragma unroll 13
	for(size_t i = 0; i < 13; ++i) {
		hmc_complex a1p;
		a1p.re = a0.re + a * a2.re;
		a1p.im = a0.im + a * a2.im;
		a0.re = fac - b * a2.im;
		a0.im =     + b * a2.re;
		a2.re = a1.re;
		a2.im = a1.im;
		a1.re = a1p.re;
		a1.im = a1p.im;
		fac *= r;
		--r;
	}

	// vr = a0 + a1*v + a2*v2
	Matrixsu3 vr;
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
	vr.e20.re =         a1.re * v.e20.re - a1.im * v.e20.im + a2.re * v2.e20.re - a2.im * v2.e20.im;
	vr.e20.im =         a1.re * v.e20.im + a1.im * v.e20.re + a2.re * v2.e20.im + a2.im * v2.e20.re;
	vr.e21.re =         a1.re * v.e21.re - a1.im * v.e21.im + a2.re * v2.e21.re - a2.im * v2.e21.im;
	vr.e21.im =         a1.re * v.e21.im + a1.im * v.e21.re + a2.re * v2.e21.im + a2.im * v2.e21.re;
	vr.e22.re = a0.re + a1.re * v.e22.re - a1.im * v.e22.im + a2.re * v2.e22.re - a2.im * v2.e22.im;
	vr.e22.im = a0.im + a1.re * v.e22.im + a1.im * v.e22.re + a2.re * v2.e22.im + a2.im * v2.e22.re;
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
	out.e20.re *= factor;
	out.e20.im *= factor;
	out.e21.re *= factor;
	out.e21.im *= factor;
	out.e22.re *= factor;
	out.e22.im *= factor;

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
	out.e20 = complexmult(in.e20, factor);
	out.e21 = complexmult(in.e21, factor);
	out.e22 = complexmult(in.e22, factor);
	return out;
}

