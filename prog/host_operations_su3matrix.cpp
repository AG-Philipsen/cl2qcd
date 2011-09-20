#include "host_operations_su3matrix.h"


//operations that contain explicit SU(3) indices!!!

hmc_complex det_su3matrix(hmc_su3matrix * U)
{
	hmc_complex det, det1, det2, det3, det4, det5, det6, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	det.re = 0;
	det.im = 0;
	tmp1 = complexmult( &(*U)[1][1], &(*U)[2][2] );
	det1 = complexmult( &(*U)[0][0] , &tmp1);
	tmp2 = complexmult( &(*U)[1][2], &(*U)[2][0] );
	det2 = complexmult( &(*U)[0][1] , &tmp2);
	tmp3 = complexmult( &(*U)[1][0], &(*U)[2][1] );
	det3 = complexmult( &(*U)[0][2] , &tmp3);
	tmp4 = complexmult( &(*U)[1][1], &(*U)[2][0] );
	det4 = complexmult( &(*U)[0][2] , &tmp4);
	tmp5 = complexmult( &(*U)[1][0], &(*U)[2][2] );
	det5 = complexmult( &(*U)[0][1] , &tmp5);
	tmp6 = complexmult( &(*U)[1][2], &(*U)[2][1] );
	det6 = complexmult( &(*U)[0][0] , &tmp6);

	det.re = det1.re + det2.re + det3.re - det4.re - det5.re - det6.re;
	det.im = det1.im + det2.im + det3.im - det4.im - det5.im - det6.im;

	return det;
}

/** @todo memcpy ... */
void copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*out)[a][b] = (*in)[a][b];
		}
	}
}

/** @todo memcpy ... */
void copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*out)[a][b] = (*in)[a][b];
		}
	}
}

/** @todo memset ... */
void zero_su3matrix(hmc_su3matrix * u)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*u)[a][b].re = 0;
			(*u)[a][b].im = 0;
		}
	}
}

/** @todo memset ... */
void zero_staplematrix(hmc_staplematrix * u)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*u)[a][b].re = 0;
			(*u)[a][b].im = 0;
		}
	}
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

void unit_su3matrix(hmc_su3matrix * u)
{
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			if(a != b) {
				(*u)[a][b].re = 0;
			} else {
				(*u)[a][b].re = hmc_one_f;
			}
			(*u)[a][b].im = 0;
		}
	}
}

void random_su3matrix(hmc_su3matrix *)
{
	throw Print_Error_Message("random su3matrix needs to be implemented...");
}

void multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q)
{
	for(int i = 0; i < NC; i++) {
		for(int k = 0; k < NC; k++) {
			(*out)[i][k].re = 0;
			(*out)[i][k].im = 0;
			for(int j = 0; j < NC; j++) {
				hmc_complex tmp = complexmult(&(*p)[i][j], &(*q)[j][k]);
				complexaccumulate(&(*out)[i][k], &tmp);
			}
		}
	}
}

void multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q)
{
	for(int i = 0; i < NC; i++) {
		for(int k = 0; k < NC; k++) {
			(*out)[i][k].re = 0;
			(*out)[i][k].im = 0;
			for(int j = 0; j < NC; j++) {
				hmc_complex tmp = complexmult(&(*p)[i][j], &(*q)[j][k]);
				complexaccumulate(&(*out)[i][k], &tmp);
			}
		}
	}
}

void accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q)
{
	for(int i = 0; i < NC; i++) {
		for(int k = 0; k < NC; k++) {
			complexaccumulate(&(*p)[i][k], &(*q)[i][k]);
		}
	}
}

void accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *mult)
{
	hmc_su3matrix tmp;
	multiply_su3matrices(&tmp, acc, mult);
	copy_su3matrix(acc, &tmp);
}

void adjoin_su3matrix(hmc_su3matrix * mat)
{
	hmc_su3matrix tmp;
	copy_su3matrix(&tmp, mat);
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			(*mat)[a][b] = complexconj(&(tmp[b][a]));
		}
	}
}

hmc_complex trace_su3matrix(hmc_su3matrix * mat)
{
	hmc_complex trace;
	trace.re = 0;
	trace.im = 0;
	for(int a = 0; a < NC; a++) complexaccumulate(&trace, &((*mat)[a][a]));;
	return trace;
}

void gaugefield_apply_bc(hmc_su3matrix * in, hmc_float theta)
{
	hmc_float tmp1, tmp2;
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			tmp1 = ((*in)[a][b]).re;
			tmp2 = ((*in)[a][b]).im;
			((*in)[a][b]).re = cos(theta) * tmp1 - sin(theta) * tmp2;
			((*in)[a][b]).im = sin(theta) * tmp1 + cos(theta) * tmp2;
		}
	}
}

// replace link in with e^(mu.re)*(cos(mu.im) + i*sin(mu.im))
void gaugefield_apply_chem_pot(hmc_su3matrix * u, hmc_su3matrix * udagger, hmc_float chem_pot_re, hmc_float chem_pot_im)
{
	hmc_float tmp1, tmp2;
	for(int a = 0; a < NC; a++) {
		for(int b = 0; b < NC; b++) {
			tmp1 = ((*u)[a][b]).re;
			tmp2 = ((*u)[a][b]).im;
			((*u)[a][b]).re = exp(chem_pot_re) * ( cos(chem_pot_im) * tmp1 - sin(chem_pot_im) * tmp2 );
			((*u)[a][b]).im = exp(chem_pot_re) * ( sin(chem_pot_im) * tmp1 + cos(chem_pot_im) * tmp2 );
			tmp1 = ((*udagger)[a][b]).re;
			tmp2 = ((*udagger)[a][b]).im;
			((*udagger)[a][b]).re = exp(-chem_pot_re) * ( cos(chem_pot_im) * tmp1 + sin(chem_pot_im) * tmp2 );
			((*udagger)[a][b]).im = exp(-chem_pot_re) * ( -sin(chem_pot_im) * tmp1 + cos(chem_pot_im) * tmp2 );
		}
	}
}

