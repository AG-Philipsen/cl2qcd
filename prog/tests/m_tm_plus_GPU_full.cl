/** @file
 * Inclusion and definition of types and definitions required in the Device code.
 */

//opencl_header.cl

#pragma OPENCL EXTENSION cl_amd_printf : enable

#ifdef _USEDOUBLEPREC_
#ifdef _DEVICE_DOUBLE_EXTENSION_AMD_
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
#ifdef _DEVICE_DOUBLE_EXTENSION_KHR_
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
#endif


#include "globaldefs.h" //NDIM, NSPIN, NC
#include "types.h"

//!!CP: why is this here?
// #define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

//for hmc_ocl_su3matrix
#ifdef _RECONSTRUCT_TWELVE_
#define SU3SIZE NC*(NC-1)
#define STAPLEMATRIXSIZE NC*NC
#else
#define SU3SIZE NC*NC
#define STAPLEMATRIXSIZE NC*NC
#endif
/** @file
 * Device code for lattice geometry handling
 */

//opencl_geometry.cl

int inline get_global_pos(int spacepos, int t)
{
return spacepos + VOLSPACE * t;
}

int inline get_global_link_pos(int mu, int spacepos, int t)
{
return mu + NDIM *get_global_pos(spacepos, t);
}

//dispensable
//site = pos + VOLSPACE*t = x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
/*
int inline ocl_gaugefield_element(int mu, int spacepos, int t)
{
return get_global_link_pos(mu, spacepos, t);
}
*/

//old version, this can be deleted soon:
int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t)
{
#ifdef _RECONSTRUCT_TWELVE_
       return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
#else
	return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*NC*NDIM*spacepos+2*NC*NC*NDIM*VOLSPACE*t;
#endif
}

int inline ocl_su3matrix_element(int a, int b)
{
#ifdef _RECONSTRUCT_TWELVE_
       return a + (NC-1)*b;
#else
	return a + NC*b;
#endif
}


//site = pos + VOLSPACE*t =  x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
// int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t)
// {
// #ifdef _RECONSTRUCT_TWELVE_
// 	return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
// #else
// 	return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*NC*NDIM*spacepos+2*NC*NC*NDIM*VOLSPACE*t;
// #endif
// }

//dispensable
/*
int inline ocl_su3matrix_element(int a, int b)
{
#ifdef _RECONSTRUCT_TWELVE_
	return a + (NC-1)*b;
#else
	return a + NC*b;
#endif
}
*/

//it is assumed that idx iterates only over half the number of sites
void inline get_even_site(int idx, int * out_space, int * out_t)
{
	int x,y,z,t;
	x = idx;
	t = (int)(idx/(VOLSPACE/2));
	x -= t*VOLSPACE/2;
	z = (int)(x/(NSPACE*NSPACE/2));
	x -= z*NSPACE*NSPACE/2;
	y = (int)(x/NSPACE);
	x -= y*NSPACE;
	(*out_space) =  (int)((z+t)%2)*(1 + 2*x - (int) (2*x/NSPACE)) + (int)((t+z+1)%2)*(2*x + (int) (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
	(*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites
void inline get_odd_site(int idx, int * out_space, int * out_t)
{
	int x,y,z,t;
	x = idx;
	t = (int)(idx/(VOLSPACE/2));
	x -= t*VOLSPACE/2;
	z = (int)(x/(NSPACE*NSPACE/2));
	x -= z*NSPACE*NSPACE/2;
	y = (int)(x/NSPACE);
	x -= y*NSPACE;
	(*out_space) =  (int)((z+t+1)%2)*(1 + 2*x - (int) (2*x/NSPACE)) + (int)((t+z)%2)*(2*x + (int) (2*x/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*z;
	(*out_t) = t;
}

int get_nspace(const int* coord)
{
	int n = coord[1] +  NSPACE*coord[2] + NSPACE*NSPACE*coord[3];
	return n;
}

int get_spacecoord(const int nspace, const int dir)
{
	int res = convert_int(nspace/(NSPACE*NSPACE));
	if(dir==3) return res;
	int acc = res;
	res = convert_int(nspace/(NSPACE)) - NSPACE*acc;
	if(dir==2) return res;
	acc = NSPACE*acc + res;
	res = nspace - NSPACE*acc;
	return res;
}

void get_allspacecoord(const int nspace, int coord[NDIM])
{
	int res = convert_int(nspace/(NSPACE*NSPACE));
	coord[3] = res;
	int acc = res;
	res = convert_int(nspace/(NSPACE)) - NSPACE*acc;
	coord[2] = res;
	acc = NSPACE*acc + res;
	res = nspace - NSPACE*acc;
	coord[1] = res;
}

int get_neighbor(const int nspace,const int dir)
{
	int coord[NDIM];
	get_allspacecoord(nspace, coord);
	coord[dir] = (coord[dir] + 1)%NSPACE;
	return get_nspace(coord);
}

int get_lower_neighbor(const int nspace, int const dir)
{
	int coord[NDIM];
	get_allspacecoord(nspace, coord);
	coord[dir] = (coord[dir] - 1 + NSPACE)%NSPACE;
	return get_nspace(coord);
}

int spinor_color(int spinor_element)
{
	return (int)(spinor_element/NSPIN);
}

int spinor_spin(int spinor_element,int color)
{
	return spinor_element - NSPIN*color;
}

int spinor_element(int alpha, int color)
{
	return alpha + NSPIN*color;
}

int get_n_eoprec(int spacepos, int timepos)
{
	return (int)((get_global_pos(spacepos, timepos))/2);
}

int eoprec_spinor_field_element(int alpha, int color, int n_eoprec)
{
	return alpha + NSPIN*color + NSPIN*NC*n_eoprec;
}

int spinor_field_element(int alpha, int color, int nspace, int t)
{
	return alpha + NSPIN*color + NSPIN*NC*(get_global_pos(nspace, t));
}
/** @file
 * Device code implementing complex numbers
 */

//opencl_operations_complex.cl

hmc_complex complexconj(hmc_complex in)
{
	in.im = -(in.im);
	return in;
}

hmc_complex complexmult(const hmc_complex a, const hmc_complex b)
{
	hmc_complex res;
	res.re = a.re * b.re - a.im * b.im;
	res.im = a.im * b.re + a.re * b.im;
	return res;
}

hmc_complex complexadd(hmc_complex a, hmc_complex b)
{
	hmc_complex res;
	res.re = a.re + b.re;
	res.im = a.im + b.im;
	return res;
}

hmc_complex complexsubtract(hmc_complex a, hmc_complex b)
{
	hmc_complex res;
	res.re = a.re - b.re;
	res.im = a.im - b.im;
	return res;
}

hmc_complex complexdivide(hmc_complex numerator, hmc_complex denominator)
{
	hmc_float norm = denominator.re * denominator.re + denominator.im * denominator.im;
	hmc_complex res;
	res.re = (numerator.re * denominator.re + numerator.im * denominator.im ) / norm;
	res.im = (numerator.im * denominator.re - numerator.re * denominator.im ) / norm;
	return res;
}

//this is deprecated...
// void gaussianComplexVector(__global hmc_complex * vector, int length, hmc_float sigma, __global hmc_ocl_ran * rnd){
// 	// SL: this fills real and imaginary part of a vector of "length" complex numbers
// 	//     with components drawn with a Gaussian distribution and variance sigma
// 	for(int idx=0;idx<length;idx++){
// 		gaussianNormalPair(&vector[idx].re,&vector[idx].im);
// 		vector[idx].re*=sigma;
// 		vector[idx].im*=sigma;
// 	}
// 	return HMC_SUCCESS;
// 	// SL: not yet tested
// }

/** @file
 * Device code implementing SU(3)matrices with and without reconstruct 12
 */
//operations_matrix_su3.cl

#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(const Matrixsu3 p, const int ncomp)
{
	hmc_complex out;
	hmc_complex tmp;
	
	switch (ncomp){

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

void print_matrixsu3(Matrixsu3 in){
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
		out.e00.im =-p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
		           + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;
	      
		out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
		           + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
		out.e01.im =-p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
		           + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;
	       
		out.e02.re = p.e00.re * q20_re   + p.e01.re * q21_re   + p.e02.re * q22_re
		           + p.e00.im * q20_im   + p.e01.im * q21_im   + p.e02.im * q22_im;
		out.e02.im =-p.e00.re * q20_im   - p.e01.re * q21_im   - p.e02.re * q22_im
		           + p.e00.im * q20_re   + p.e01.im * q21_re   + p.e02.im * q22_re;
	       
		out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
		           + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
		out.e10.im =-p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
		           + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;
	       
		out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
		           + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
		out.e11.im =-p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
		           + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;
	       
		out.e12.re = p.e10.re * q20_re   + p.e11.re * q21_re   + p.e12.re * q22_re
		           + p.e10.im * q20_im   + p.e11.im * q21_im   + p.e12.im * q22_im;
		out.e12.im =-p.e10.re * q20_im   - p.e11.re * q21_im   - p.e12.re * q22_im
		           + p.e10.im * q20_re   + p.e11.im * q21_re   + p.e12.im * q22_re;	       
#else
    
		out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
		           + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
		out.e00.im =-p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
		           + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

		out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
		           + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
		out.e01.im =-p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
		           + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;
    
		out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
		           + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
		out.e02.im =-p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
		           + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;
	       
		out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
		           + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
		out.e10.im =-p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
		           + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;
	       
		out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
		           + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
		out.e11.im =-p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
		           + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;
	       
		out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
		           + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
		out.e12.im =-p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
		           + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;	       

		out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
		           + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
		out.e20.im =-p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
		           + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;
	       
		out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
		           + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
		out.e21.im =-p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
		           + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;
	       
		out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
		           + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
		out.e22.im =-p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
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
Matrixsu3 build_su3_from_ae(ae in){
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
Matrixsu3 build_su3_from_ae_times_real(ae in, hmc_float eps){
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
  Matrixsu3 v,v2,vr;
  hmc_float fac,r;
  hmc_float a,b;
  hmc_complex a0,a1,a2,a1p;

  //original in tmlqcd: _make_su3(v,p);
	//Here, the stepsize factor is multiplied in right away!!
	//Also, a factor of 0.5 is taken out to fit the different trace-definition from tmqlcd
	hmc_float halfeps = epsilon;// *F_1_2;
	v = build_su3_from_ae_times_real(inn, halfeps);
	
  // calculates v^2 
  v2 = multiply_matrixsu3(v,v);
#ifdef _RECONSTRUCT_TWELVE_
	hmc_complex v_e20 = reconstruct_su3(v, 0);
	hmc_complex v_e21 = reconstruct_su3(v, 1);
	hmc_complex v_e22 = reconstruct_su3(v, 2);
	
	hmc_complex v2_e20 = reconstruct_su3(v2, 0);
	hmc_complex v2_e21 = reconstruct_su3(v2, 1);
	hmc_complex v2_e22 = reconstruct_su3(v2, 2);
	
	a=0.5*(v2.e00.re+v2.e11.re+v2_e22.re);
	
	b = 0.33333333333333333*
    (v.e00.re*v2.e00.im+v.e00.im*v2.e00.re
     +v.e01.re*v2.e10.im+v.e01.im*v2.e10.re
     +v.e02.re*v2_e20.im+v.e02.im*v2_e20.re
     +v.e10.re*v2.e01.im+v.e10.im*v2.e01.re
     +v.e11.re*v2.e11.im+v.e11.im*v2.e11.re
     +v.e12.re*v2_e21.im+v.e12.im*v2_e21.re
     +v_e20.re*v2.e02.im+v_e20.im*v2.e02.re
     +v_e21.re*v2.e12.im+v_e21.im*v2.e12.re
     +v_e22.re*v2_e22.im+v_e22.im*v2_e22.re  );
	
#else
  a=0.5*(v2.e00.re+v2.e11.re+v2.e22.re);
  // 1/3 imaginary part of tr v*v2 
  b = 0.33333333333333333*
    (v.e00.re*v2.e00.im+v.e00.im*v2.e00.re
     +v.e01.re*v2.e10.im+v.e01.im*v2.e10.re
     +v.e02.re*v2.e20.im+v.e02.im*v2.e20.re
     +v.e10.re*v2.e01.im+v.e10.im*v2.e01.re
     +v.e11.re*v2.e11.im+v.e11.im*v2.e11.re
     +v.e12.re*v2.e21.im+v.e12.im*v2.e21.re
     +v.e20.re*v2.e02.im+v.e20.im*v2.e02.re
     +v.e21.re*v2.e12.im+v.e21.im*v2.e12.re
     +v.e22.re*v2.e22.im+v.e22.im*v2.e22.re  );
#endif
  a0.re=0.16059043836821615e-9;    //  1/13! 
  a0.im=0.0;
  a1.re=0.11470745597729725e-10;   //  1/14! 
  a1.im=0.0;
  a2.re=0.76471637318198165e-12;   //  1/15! 
  a2.im=0.0;
  fac=0.20876756987868099e-8;      //  1/12! 
  r=12.0;
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
  vr.e00.re = a0.re + a1.re*v.e00.re - a1.im*v.e00.im + a2.re*v2.e00.re - a2.im*v2.e00.im;
  vr.e00.im = a0.im + a1.re*v.e00.im + a1.im*v.e00.re + a2.re*v2.e00.im + a2.im*v2.e00.re;
  vr.e01.re =         a1.re*v.e01.re - a1.im*v.e01.im + a2.re*v2.e01.re - a2.im*v2.e01.im;
  vr.e01.im =         a1.re*v.e01.im + a1.im*v.e01.re + a2.re*v2.e01.im + a2.im*v2.e01.re;
  vr.e02.re =         a1.re*v.e02.re - a1.im*v.e02.im + a2.re*v2.e02.re - a2.im*v2.e02.im;
  vr.e02.im =         a1.re*v.e02.im + a1.im*v.e02.re + a2.re*v2.e02.im + a2.im*v2.e02.re;
  vr.e10.re =         a1.re*v.e10.re - a1.im*v.e10.im + a2.re*v2.e10.re - a2.im*v2.e10.im;
  vr.e10.im =         a1.re*v.e10.im + a1.im*v.e10.re + a2.re*v2.e10.im + a2.im*v2.e10.re;
  vr.e11.re = a0.re + a1.re*v.e11.re - a1.im*v.e11.im + a2.re*v2.e11.re - a2.im*v2.e11.im;
  vr.e11.im = a0.im + a1.re*v.e11.im + a1.im*v.e11.re + a2.re*v2.e11.im + a2.im*v2.e11.re;
  vr.e12.re =         a1.re*v.e12.re - a1.im*v.e12.im + a2.re*v2.e12.re - a2.im*v2.e12.im;
  vr.e12.im =         a1.re*v.e12.im + a1.im*v.e12.re + a2.re*v2.e12.im + a2.im*v2.e12.re;
#ifndef _RECONSTRUCT_TWELVE_
	vr.e20.re =         a1.re*v.e20.re - a1.im*v.e20.im + a2.re*v2.e20.re - a2.im*v2.e20.im;
  vr.e20.im =         a1.re*v.e20.im + a1.im*v.e20.re + a2.re*v2.e20.im + a2.im*v2.e20.re;
  vr.e21.re =         a1.re*v.e21.re - a1.im*v.e21.im + a2.re*v2.e21.re - a2.im*v2.e21.im;
  vr.e21.im =         a1.re*v.e21.im + a1.im*v.e21.re + a2.re*v2.e21.im + a2.im*v2.e21.re;
  vr.e22.re = a0.re + a1.re*v.e22.re - a1.im*v.e22.im + a2.re*v2.e22.re - a2.im*v2.e22.im;
  vr.e22.im = a0.im + a1.re*v.e22.im + a1.im*v.e22.re + a2.re*v2.e22.im + a2.im*v2.e22.re;
#endif
  return vr;
	
}

//scale a su3 matrix by a real factor
Matrixsu3 multiply_matrixsu3_by_real (Matrixsu3 in, hmc_float factor){
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

Matrixsu3 multiply_matrixsu3_by_complex (Matrixsu3 in, hmc_complex factor){
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

/** @file
 * Device code implementing 3x3 matrices
 */

//operations_matrix.cl


void print_matrix3x3(Matrix3x3 in){
	printf("(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n(%f,%f) (%f,%f) (%f,%f)\n", 
					in.e00.re, in.e00.im, in.e01.re, in.e01.im, in.e02.re, in.e02.im, 
					in.e10.re, in.e10.im, in.e11.re, in.e11.im, in.e12.re, in.e12.im, 
					in.e20.re, in.e20.im, in.e21.re, in.e21.im, in.e22.re, in.e22.im);
	printf("\n");
}

Matrix3x3 zero_matrix3x3 (){
    Matrix3x3 out;
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

Matrix3x3 identity_matrix3x3 (){
    Matrix3x3 out;
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

Matrix3x3 multiply_matrix3x3_by_real (Matrix3x3 in, hmc_float factor){
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

Matrix3x3 multiply_matrix3x3_by_complex (Matrix3x3 in, hmc_complex factor){
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

Matrix3x3 multiply_matrix3x3 (const Matrix3x3 p, const Matrix3x3 q)
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

Matrix3x3 multiply_matrix3x3_dagger(const Matrix3x3 p, const Matrix3x3 q)
{
    Matrix3x3 out;
		
		out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
		           + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
		out.e00.im =-p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
		           + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

		out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
		           + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
		out.e01.im =-p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
		           + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;
    
		out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
		           + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
		out.e02.im =-p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
		           + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;
	       
		out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
		           + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
		out.e10.im =-p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
		           + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;
	       
		out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
		           + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
		out.e11.im =-p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
		           + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;
	       
		out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
		           + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
		out.e12.im =-p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
		           + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;	       

		out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
		           + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
		out.e20.im =-p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
		           + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;
	       
		out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
		           + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
		out.e21.im =-p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
		           + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;
	       
		out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
		           + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
		out.e22.im =-p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
		           + p.e20.im * q.e20.re + p.e21.im * q.e21.re + p.e22.im * q.e22.re;
    return out;
}



Matrix3x3 add_matrix3x3 ( const Matrix3x3 p, const Matrix3x3 q)
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

Matrix3x3 subtract_matrix3x3 ( const Matrix3x3 p, const Matrix3x3 q)
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

Matrix3x3 copy_matrix3x3(const Matrix3x3 p)
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

hmc_complex trace_matrix3x3 (const Matrix3x3 p)
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

Matrix3x3 adjoint_matrix3x3 (const Matrix3x3 p){

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

Matrix3x3 matrix_su3to3x3 (const Matrixsu3 p){

	Matrix3x3 out;
	out.e00.re = p.e00.re;
	out.e00.im = p.e00.im;
	out.e01.re = p.e01.re;
	out.e01.im = p.e01.im;
	out.e02.re = p.e02.re;
	out.e02.im = p.e02.im;

	out.e10.re = p.e10.re;
	out.e10.im = p.e10.im;
	out.e11.re = p.e11.re;
	out.e11.im = p.e11.im;
	out.e12.re = p.e12.re;
	out.e12.im = p.e12.im;

#ifdef _RECONSTRUCT_TWELVE_
	out.e20.re = reconstruct_su3(p, 0).re;
	out.e20.im = reconstruct_su3(p, 0).im;
	out.e21.re = reconstruct_su3(p, 1).re;
	out.e21.im = reconstruct_su3(p, 1).im;
	out.e22.re = reconstruct_su3(p, 2).re;
	out.e22.im = reconstruct_su3(p, 2).im;
#else
	out.e20.re = p.e20.re;
	out.e20.im = p.e20.im;
	out.e21.re = p.e21.re;
	out.e21.im = p.e21.im;
	out.e22.re = p.e22.re;
	out.e22.im = p.e22.im;
#endif

	return out;
}

hmc_float absoluteDifference_matrix3x3(Matrix3x3 mat1, Matrix3x3 mat2)
{
	hmc_float result=0.0;
	
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

Matrixsu3 project_anti_herm(Matrix3x3 in) {
  Matrixsu3 tmp;
	hmc_float tr_omega ; 

	//in is not a su3 matrix, so reconstruct does not have to be applied here!!
	tr_omega = (in.e00.im + in.e11.im + in.e22.im) / 3.0 ; 
	
	tmp.e00.re = 0.0 ;
  tmp.e00.im = -in.e00.im + tr_omega ;

  tmp.e11.re = 0.0 ;
  tmp.e11.im = -in.e11.im + tr_omega ;
	
#ifndef _RECONSTRUCT_TWELVE_
  tmp.e22.re = 0.0 ;
  tmp.e22.im = -in.e22.im + tr_omega ;
#endif

	// is this all right? check with original again!!
  tmp.e01.re = in.e01.re - in.e10.re ;  tmp.e01.re /= -2.0 ; 
  tmp.e01.im = in.e01.im + in.e10.im ;  tmp.e01.im /= -2.0 ; 
  tmp.e10.re  = -tmp.e01.re ; 
  tmp.e10.im  =  tmp.e01.im ; 

  tmp.e02.re = in.e01.re - in.e20.re ;  tmp.e02.re /= -2.0 ; 
  tmp.e02.im = in.e02.im + in.e20.im ;  tmp.e02.im /= -2.0 ; 
#ifndef _RECONSTRUCT_TWELVE_
	tmp.e20.re  = -tmp.e02.re ; 
  tmp.e20.im  =  tmp.e02.im ; 
#endif
  tmp.e12.re = in.e12.re - in.e21.re ;  tmp.e12.re /= -2.0 ; 
  tmp.e12.im = in.e12.im + in.e21.im ;  tmp.e12.im /= -2.0 ; 
#ifndef _RECONSTRUCT_TWELVE_
	tmp.e21.re  = -tmp.e12.re ; 
  tmp.e21.im  =  tmp.e12.im ; 
#endif
  // flip sign, this is just taken from the tmlqcd analogue
  tmp = multiply_matrixsu3_by_real(tmp, -1.0) ; 
	
	return tmp;

}
/** @file
 * Device code for operations on the fermion matrix
 */

//operations_gaugefield.cl


// hmc_complex global_trace_su3(__global hmc_ocl_gaugefield * field, int mu)
// {
//  hmc_complex sum;
//  sum.re = 0;
//  sum.im = 0;
//  for(int t=0; t<NTIME; t++) {
//    for(int n=0; n<VOLSPACE; n++) {
//      hmc_ocl_su3matrix tmp[SU3SIZE];
//      get_su3matrix(tmp, field, n, t, mu);
//      sum.re += trace_su3matrix(tmp).re;
//      sum.im += trace_su3matrix(tmp).im;
//    }
//  }
//  return sum;
// }

Matrixsu3 project_su3(const Matrixsu3 U)
{

	Matrixsu3 out;

	//Extract initial vectors
	hmc_complex a[NC];
	hmc_complex b[NC];

#ifdef _RECONSTRUCT_TWELVE_
	a[0] = U.e00;
	a[1] = U.e01;
	a[2] = U.e02;
	b[0] = U.e10;
	b[1] = U.e11;
	b[2] = U.e12;
#else
	hmc_complex c[NC];
	a[0] = U.e00;
	a[1] = U.e01;
	a[2] = U.e02;
	b[0] = U.e10;
	b[1] = U.e11;
	b[2] = U.e12;
	c[0] = U.e20;
	c[1] = U.e21;
	c[2] = U.e22;
#endif

	//New SU3-Matrix
	//first vector
	//norm
	hmc_float norm = 0.;
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp = complexconj(a[i]);
		tmp = complexmult (a[i], tmp);
		norm += tmp.re;
	}
	norm = 1. / sqrt(norm);
	//rescale
	for (int i = 0; i < NC; i++) {
		a[i].re *= norm;
		a[i].im *= norm;
	}

	//second vector
	//orthogonal vector
	hmc_complex factor;
	factor.re = 0.0;
	factor.im = 0.0;
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp;
		tmp = complexconj (b[i]);
		tmp = complexmult (a[i], tmp);
		factor = complexadd (factor, tmp);
	}
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp;
		tmp = complexmult(factor, a[i]);
		b[i] = complexsubtract(b[i], tmp);
	}

//norm
	norm = 0.;
	for (int i = 0; i < NC; i++) {
		hmc_complex tmp;
		tmp = complexconj(b[i]);
		tmp = complexmult (b[i], tmp);
		norm +=  tmp.re;
	}
	norm = 1. / sqrt(norm);
	//rescale
	for  (int i = 0; i < NC; i++) {
		b[i].re *= norm;
		b[i].im *= norm;
	}

#ifndef _RECONSTRUCT_TWELVE_
	//third vector
	//orthogonal vector
	hmc_complex tmp;
	hmc_complex tmp2;
	tmp = complexmult(a[1], b[2]);
	tmp = complexconj(tmp);
	tmp2 = complexmult(a[2], b[1]);
	tmp2 = complexconj(tmp2);
	c[0] = complexsubtract(tmp, tmp2);
	tmp = complexmult(a[2], b[0]);
	tmp = complexconj(tmp);
	tmp2 = complexmult(a[0], b[2]);
	tmp2 = complexconj(tmp2);
	c[1] = complexsubtract(tmp, tmp2);
	tmp = complexmult(a[0], b[1]);
	tmp = complexconj(tmp);
	tmp2 = complexmult(a[1], b[0]);
	tmp2 = complexconj(tmp2);
	c[2] = complexsubtract(tmp, tmp2);
	//Set new values to matrix
	out.e20 = c[0];
	out.e21 = c[1];
	out.e22 = c[2];
#endif

	//Set new values to matrix
	out.e00 = a[0];
	out.e10 = b[0];
	out.e01 = a[1];
	out.e11 = b[1];
	out.e02 = a[2];
	out.e12 = b[2];

	return out;
}

Matrixsu2 reduction (const Matrix3x3 src, const int rand)
{
	Matrixsu2 out;
	if(rand == 1) {
		out.e00 = src.e00;
		out.e01 = src.e01;
		out.e10 = src.e10;
		out.e11 = src.e11;
	} else if (rand == 2) {
		out.e00 = src.e11;
		out.e01 = src.e12;
		out.e10 = src.e21;
		out.e11 = src.e22;
	} else if (rand == 3) {
		out.e00 = src.e00;
		out.e01 = src.e02;
		out.e10 = src.e20;
		out.e11 = src.e22;
	}
	return out;
}

//CP: I leave both version in here...
Matrixsu3 extend (const int random, Matrixsu2 src)
{
	Matrixsu3 out;

#ifdef _RECONSTRUCT_TWELVE_
	switch(random) {
		case 1:
			out.e00 = src.e00;
			out.e01 = src.e01;
			out.e02 = hmc_complex_zero;
			out.e10 = src.e10;
			out.e11 = src.e11;
			out.e12 = hmc_complex_zero;
			return out;
		case 2:
			out.e00 = hmc_complex_one;
			out.e01 = hmc_complex_zero;
			out.e02 = hmc_complex_zero;
			out.e10 = hmc_complex_zero;
			out.e11 = src.e00;
			out.e12 = src.e01;
			return out;
		case 3:
			out.e00 = src.e00;
			out.e01 = hmc_complex_zero;
			out.e02 = src.e01;
			out.e10 = hmc_complex_zero;
			out.e11 = hmc_complex_one;
			out.e12 = hmc_complex_zero;
			return out;
	}
#else

	switch(random) {
		case 1:
			out.e00 = src.e00;
			out.e01 = src.e01;
			out.e02 = hmc_complex_zero;
			out.e10 = src.e10;
			out.e11 = src.e11;
			out.e12 = hmc_complex_zero;
			out.e20 = hmc_complex_zero;
			out.e21 = hmc_complex_zero;
			out.e22 = hmc_complex_one;
			return out;
		case 2:
			out.e00 = hmc_complex_one;
			out.e01 = hmc_complex_zero;
			out.e02 = hmc_complex_zero;
			out.e10 = hmc_complex_zero;
			out.e11 = src.e00;
			out.e12 = src.e01;
			out.e20 = hmc_complex_zero;
			out.e21 = src.e10;
			out.e22 = src.e11;
			return out;
		case 3:
			out.e00 = src.e00;
			out.e01 = hmc_complex_zero;
			out.e02 = src.e01;
			out.e10 = hmc_complex_zero;
			out.e11 = hmc_complex_one;
			out.e12 = hmc_complex_zero;
			out.e20 = src.e10;
			out.e21 = hmc_complex_zero;
			out.e22 = src.e11;
			return out;
	}
#endif
	return;
}

// void gaugefield_apply_bc(__private hmc_ocl_su3matrix * in, hmc_float theta)
// {
//  hmc_float tmp1,tmp2;
// #ifdef _RECONSTRUCT_TWELVE_
//  for(int a=0; a<NC; a++) {
//    for(int b=0; b<NC-1; b++) {
// //       for(int n=0; n<NC*(NC-1); n++) {
//      tmp1 = ((in)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((in)[ocl_su3matrix_element(a,b)]).im;
//      ((in)[ocl_su3matrix_element(a,b)]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
//      ((in)[ocl_su3matrix_element(a,b)]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
//    }
//  }
// #else
//  for(int a=0; a<NC; a++) {
//    for(int b=0; b<NC; b++) {
//      tmp1 = ((in)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((in)[ocl_su3matrix_element(a,b)]).im;
//      ((in)[ocl_su3matrix_element(a,b)]).re = cos(theta)*tmp1 - sin(theta)*tmp2;
//      ((in)[ocl_su3matrix_element(a,b)]).im = sin(theta)*tmp1 + cos(theta)*tmp2;
//    }
//  }
// #endif
//  return;
// }

//CP: the following two chem-pot functions are one in the host-code
// void gaugefield_apply_chem_pot_real(__private hmc_ocl_su3matrix * u, __private hmc_ocl_su3matrix * udagger, hmc_float chem_pot_re)
// {
//  hmc_float tmp1,tmp2;
// #ifdef _RECONSTRUCT_TWELVE_
//  for(int a=0; a<NC; a++) {
//    for(int b=0; b<NC-1; b++) {
// //   for(int n=0; n<NC*(NC-1); n++) {
//      tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
//      ((u)[ocl_su3matrix_element(a,b)]).re = exp(chem_pot_re)*tmp1;
//      ((u)[ocl_su3matrix_element(a,b)]).im = exp(chem_pot_re)*tmp2;
//      tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
//      ((udagger)[ocl_su3matrix_element(a,b)]).re = exp(-chem_pot_re)*tmp1;
//      ((udagger)[ocl_su3matrix_element(a,b)]).im = exp(-chem_pot_re)*tmp2;
//    }
//  }
// #else
//  for(int a=0; a<NC; a++) {
//    for(int b=0; b<NC; b++) {
//      tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
//      ((u)[ocl_su3matrix_element(a,b)]).re = exp(chem_pot_re)*tmp1;
//      ((u)[ocl_su3matrix_element(a,b)]).im = exp(chem_pot_re)*tmp2;
//      tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
//      ((udagger)[ocl_su3matrix_element(a,b)]).re = exp(-chem_pot_re)*tmp1;
//      ((udagger)[ocl_su3matrix_element(a,b)]).im = exp(-chem_pot_re)*tmp2;
//    }
//  }
// #endif
//  return;
// }

// void gaugefield_apply_chem_pot_imag(__private hmc_ocl_su3matrix * u, __private hmc_ocl_su3matrix * udagger, hmc_float chem_pot_im)
// {
//  hmc_float tmp1,tmp2;
// #ifdef _RECONSTRUCT_TWELVE_
//  for(int a=0; a<NC; a++) {
//    for(int b=0; b<NC-1; b++) {
// //   for(int n=0; n<NC*(NC-1); n++) {
//      tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
//      ((u)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
//      ((u)[ocl_su3matrix_element(a,b)]).im = ( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
//      tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
//      ((udagger)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
//      ((udagger)[ocl_su3matrix_element(a,b)]).im = ( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
//    }
//  }
// #else
//  for(int a=0; a<NC; a++) {
//    for(int b=0; b<NC; b++) {
//      tmp1 = ((u)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((u)[ocl_su3matrix_element(a,b)]).im;
//      ((u)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 - sin(chem_pot_im)*tmp2 );
//      ((u)[ocl_su3matrix_element(a,b)]).im = ( sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
//      tmp1 = ((udagger)[ocl_su3matrix_element(a,b)]).re;
//      tmp2 = ((udagger)[ocl_su3matrix_element(a,b)]).im;
//      ((udagger)[ocl_su3matrix_element(a,b)]).re = ( cos(chem_pot_im)*tmp1 + sin(chem_pot_im)*tmp2 );
//      ((udagger)[ocl_su3matrix_element(a,b)]).im = ( -sin(chem_pot_im)*tmp1 + cos(chem_pot_im)*tmp2 );
//    }
//  }
// #endif
//  return;
// }

//calculate polyakov-loop matrix at spatial site n in time-direction
Matrixsu3 local_polyakov(__global ocl_s_gaugefield * field, const int n)
{
	Matrixsu3 out;
	out = unit_matrixsu3();
	for(int t = 0; t < NTIME; t++) {
		Matrixsu3 tmp;
		tmp = get_matrixsu3(field, n, t, 0);
		out = multiply_matrixsu3 (out, tmp);
	}
	return out;
}

//CP: final version, using one third of the original scratch-registers...
//calculate plaquette-matrix at site n,t in direction mu and nu
Matrixsu3 local_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	Matrixsu3 out;
	int4 pos;
	if(mu == 0) {
		pos.x = (t + 1) % NTIME;
		pos.y = n;
	} else {
		pos.x = t;
		pos.y = get_neighbor(n, mu);
	}
	if(nu == 0) {
		pos.z = (t + 1) % NTIME;
		pos.w = n;
	} else {
		pos.z = t;
		pos.w = get_neighbor(n, nu);
	}
	out = multiply_matrixsu3 (get_matrixsu3(field, n, t, mu), get_matrixsu3(field, pos.y, pos.x, nu)      );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, pos.w, pos.z, mu) );
	out = multiply_matrixsu3_dagger(out, get_matrixsu3(field, n, t, nu) );

	return out;
}


//todo
Matrix3x3 local_Q_plaquette(__global ocl_s_gaugefield * field, const int n, const int t, const int mu, const int nu )
{
	Matrix3x3 out;
	Matrixsu3 tmp;
	int newpos;

	//1st plaq
	Matrixsu3 plaq1;
	//u_mu(x)
	plaq1 = get_matrixsu3(field, n, t, mu);
	//u_nu(x+mu)
	if(mu == 0) {
		int newt = (t + 1) % NTIME;
		tmp = get_matrixsu3(field, n, newt, nu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(n, mu), t, nu);
	plaq1 = multiply_matrixsu3 (plaq1, tmp);
	//adjoint(u_mu(x+nu))
	if(nu == 0) {
		int newt = (t + 1) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(n, nu), t, mu);
	tmp = adjoint_matrixsu3(tmp);
	plaq1 = multiply_matrixsu3 (plaq1, tmp);
	//adjoint(u_nu(x))
	tmp = get_matrixsu3(field, n, t, nu);
	tmp = adjoint_matrixsu3(tmp);
	plaq1 = multiply_matrixsu3 (plaq1, tmp);

	//2nd plaq
	Matrixsu3 plaq2;
	//U_nu(x)
	plaq2 = get_matrixsu3(field, n, t, nu);
	//adj (u_mu(x-mu+nu))
	newpos = get_lower_neighbor(n, mu);
	if (nu == 0) {
		int newt =  (t + 1) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, mu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_neighbor(n, nu), newt, mu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(newpos, nu), t, mu);

	tmp = adjoint_matrixsu3(tmp);
	plaq2 = multiply_matrixsu3 (plaq2, tmp);
	//adj (u_nu(x-mu))
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, nu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, mu), t, nu);
	tmp = adjoint_matrixsu3(tmp);
	plaq2 = multiply_matrixsu3 (plaq2, tmp);
	//u_mu(x-mu)
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, mu), t, mu);
	plaq2 = multiply_matrixsu3 (plaq2, tmp);

	//3rd plaq
	Matrixsu3 plaq3;
	//adj (u_mu(x-mu))
	if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, mu), t, mu);
	tmp = adjoint_matrixsu3(tmp);
	plaq3 = copy_matrixsu3(tmp);
	//adj (u_nu(x-mu-nu))
	newpos = get_lower_neighbor(n, mu);
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, nu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), newt, nu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(newpos, nu), t, nu);
	tmp = adjoint_matrixsu3(tmp);
	plaq3 = multiply_matrixsu3 (plaq3, tmp);
	//u_mu(x-mu-nu)
	newpos = get_lower_neighbor(n, mu);
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, mu);
	} else if (mu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), newt, mu);
	} else
		tmp = get_matrixsu3 (field, get_lower_neighbor(newpos, nu), t, mu);
	plaq3 = multiply_matrixsu3(plaq3, tmp);
	//u_nu(x-nu)
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, nu);
	} else
		//this causes for nu=1 a speicherzugriffsfehler
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), t, nu);
	plaq3 = multiply_matrixsu3 (plaq3, tmp);

	//4th plaq
	Matrixsu3 plaq4;
	//adj(u_nu(x-nu))
	tmp = adjoint_matrixsu3(tmp);
	plaq4 = copy_matrixsu3(tmp);
	//u_mu(x-nu)
	if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, n, newt, mu);
	} else
		tmp = get_matrixsu3(field, get_lower_neighbor(n, nu), t, mu);
	plaq4 = multiply_matrixsu3 (plaq4, tmp);
	//u_nu(x+mu-nu)
	newpos = get_lower_neighbor(n, nu);
	if (mu == 0) {
		int newt =  (t + 1) % NTIME;
		tmp = get_matrixsu3(field, newpos, newt, nu);
	} else if (nu == 0) {
		int newt = (t - 1 + NTIME) % NTIME;
		tmp = get_matrixsu3(field, get_neighbor(n, mu), newt, nu);
	} else
		tmp = get_matrixsu3(field, get_neighbor(newpos, mu), t, nu);
	plaq4 = multiply_matrixsu3(plaq4, tmp);
	//adj (u_mu(x))
	tmp = get_matrixsu3(field, n, t, mu);
	tmp = adjoint_matrixsu3(tmp);
	plaq4 = multiply_matrixsu3 (plaq4, tmp);

	//Sum up
	Matrix3x3 tmp3x3;
	out = matrix_su3to3x3 (plaq1);
	tmp3x3 = matrix_su3to3x3 (plaq2);
	out = add_matrix3x3 (out, tmp3x3);
	tmp3x3 = matrix_su3to3x3 (plaq3);
	out = add_matrix3x3 (out, tmp3x3);
	tmp3x3 = matrix_su3to3x3 (plaq4);
	out = add_matrix3x3 (out, tmp3x3);

	return out;
}

Matrix3x3 calc_staple(__global ocl_s_gaugefield* field, const int pos, const int t, const int mu_in)
{
	Matrixsu3 prod;
	Matrixsu3 prod2;
	Matrixsu3 tmp;
	Matrix3x3 staple;
	int nu, newpos, newt;

	staple = zero_matrix3x3();

	//iterate through the three directions other than mu
	for(int i = 1; i < NDIM; i++) {

		prod = zero_matrixsu3();

		nu = (mu_in + i) % NDIM;
		//first staple
		//u_nu(x+mu)
		if(mu_in == 0) {
			newt = (t + 1) % NTIME;
			tmp = get_matrixsu3(field, pos, newt, nu);

		} else {
			tmp = get_matrixsu3(field, get_neighbor(pos, mu_in), t, nu);
		}
		prod = copy_matrixsu3(tmp);

		//adjoint(u_mu(x+nu))
		if(nu == 0) {
			newt = (t + 1) % NTIME;
			tmp = get_matrixsu3(field, pos, newt, mu_in);
		} else {
			tmp = get_matrixsu3(field, get_neighbor(pos, nu), t, mu_in);
		}
		tmp = adjoint_matrixsu3(tmp);

		prod = multiply_matrixsu3(prod, tmp);

		//adjoint(u_nu(x))
		tmp = get_matrixsu3(field, pos, t, nu);
		tmp = adjoint_matrixsu3(tmp);
		prod = multiply_matrixsu3 (prod, tmp);

		//second staple
		//adjoint (u_nu(x+mu-nu))
		//newpos is "pos-nu" (spatial)
		newpos = get_lower_neighbor(pos, nu);
		if(mu_in == 0) {
			newt = (t + 1) % NTIME;
			tmp = get_matrixsu3(field, newpos, newt, nu);
		} else if (nu == 0) {
			newt = (t - 1 + NTIME) % NTIME;
			tmp = get_matrixsu3(field, get_neighbor(pos, mu_in), newt, nu);
		} else {
			tmp = get_matrixsu3(field, get_neighbor(newpos, mu_in), t, nu);
		}
		prod2 = adjoint_matrixsu3(tmp);
		//adjoint(u_mu(x-nu))
		if(mu_in == 0) {
			tmp = get_matrixsu3(field, newpos, t, mu_in);
		} else if (nu == 0) {
			newt = (t - 1 + NTIME) % NTIME;
			tmp = get_matrixsu3(field, pos, newt, mu_in);
		} else {
			tmp = get_matrixsu3(field, newpos, t, mu_in);
		}
		tmp = adjoint_matrixsu3(tmp);
		prod2 = multiply_matrixsu3(prod2, tmp);
		//adjoint(u_nu(x-nu))
		if(mu_in == 0) {
			tmp = get_matrixsu3(field, newpos, t, nu);
		} else if (nu == 0) {
			newt = (t - 1 + NTIME) % NTIME;
			tmp = get_matrixsu3(field, pos, newt, nu);
		} else {
			tmp = get_matrixsu3(field, newpos, t, nu);
		}
		prod2 = multiply_matrixsu3(prod2, tmp);


		Matrix3x3 dummy;
		dummy = matrix_su3to3x3 (prod);
		staple = add_matrix3x3 (staple, dummy );
		dummy = matrix_su3to3x3 (prod2);
		staple = add_matrix3x3 (staple, dummy );
	}

	return staple;
}
/** @file
 * Common fermion types used by HMC, both on host and OpenCL device.
 */

#ifndef _TYPES_FERMIONSH_
#define _TYPES_FERMIONSH_

#include "types.h"

//CP: new defs for spinors
typedef struct {
	hmc_complex e0;
	hmc_complex e1;
	hmc_complex e2;
	
} su3vec;

typedef struct {
	su3vec e0;
	su3vec e1;
	su3vec e2;
	su3vec e3;
} spinor;

typedef struct {
	su3vec e0;
	su3vec e1;
} halfspinor;

typedef spinor spinorfield;
typedef spinor spinorfield_eoprec;

#endif /* _TYPES_FERMIONSH_ */


/** @file
 * Device code for operations on SU(3) vectors
 */

hmc_float su3vec_squarenorm(su3vec in){
	return
		in.e0.re*in.e0.re + in.e0.im*in.e0.im + 
		in.e1.re*in.e1.re + in.e1.im*in.e1.im + 
		in.e2.re*in.e2.re + in.e2.im*in.e2.im;
}

hmc_complex su3vec_scalarproduct(su3vec in1, su3vec in2){
	    hmc_complex tmp;
	tmp.re = 
		in1.e0.re*in2.e0.re + in1.e0.im*in2.e0.im + 
		in1.e1.re*in2.e1.re + in1.e1.im*in2.e1.im + 
		in1.e2.re*in2.e2.re + in1.e2.im*in2.e2.im;
	tmp.im = 
		in1.e0.re*in2.e0.im - in2.e0.re*in1.e0.im + 
		in1.e1.re*in2.e1.im - in2.e1.re*in1.e1.im + 
		in1.e2.re*in2.e2.im - in2.e2.re*in1.e2.im;
	return tmp;
}

su3vec set_su3vec_zero(){
	su3vec tmp;
  (tmp).e0 = hmc_complex_zero;
	(tmp).e1 = hmc_complex_zero;
	(tmp).e2 = hmc_complex_zero;
  return tmp;
}

su3vec set_su3vec_cold(){
	su3vec tmp;
  (tmp).e0 = hmc_complex_one;
	(tmp).e1 = hmc_complex_one;
	(tmp).e2 = hmc_complex_one;
  return tmp;
}

su3vec su3vec_times_real(su3vec in, hmc_float factor){
	su3vec tmp;
	tmp.e0.re = in.e0.re*factor;
	tmp.e0.im = in.e0.im*factor;
	tmp.e1.re = in.e1.re*factor;
	tmp.e1.im = in.e1.im*factor;
	tmp.e2.re = in.e2.re*factor;
	tmp.e2.im = in.e2.im*factor;
	return tmp;
}

su3vec su3vec_times_complex(su3vec in, hmc_complex factor){
	su3vec tmp;
	tmp.e0.re = in.e0.re*factor.re - in.e0.im*factor.im;
	tmp.e0.im = in.e0.im*factor.re + in.e0.re*factor.im;
	tmp.e1.re = in.e1.re*factor.re - in.e1.im*factor.im;
	tmp.e1.im = in.e1.im*factor.re + in.e1.re*factor.im;
	tmp.e2.re = in.e2.re*factor.re - in.e2.im*factor.im;
	tmp.e2.im = in.e2.im*factor.re + in.e2.re*factor.im;
	return tmp;
}

su3vec su3vec_times_complex_conj(su3vec in, hmc_complex factor){
	su3vec tmp;
	tmp.e0.re = in.e0.re*factor.re + in.e0.im*factor.im;
	tmp.e0.im = in.e0.im*factor.re - in.e0.re*factor.im;
	tmp.e1.re = in.e1.re*factor.re + in.e1.im*factor.im;
	tmp.e1.im = in.e1.im*factor.re - in.e1.re*factor.im;
	tmp.e2.re = in.e2.re*factor.re + in.e2.im*factor.im;
	tmp.e2.im = in.e2.im*factor.re - in.e2.re*factor.im;
	return tmp;
}

su3vec su3matrix_times_su3vec(Matrixsu3 u, su3vec in){
	su3vec tmp;
	#ifdef _RECONSTRUCT_TWELVE_
	hmc_float u_e20_re = reconstruct_su3(u, 0).re;
	hmc_float u_e20_im = reconstruct_su3(u, 0).im;
	hmc_float u_e21_re = reconstruct_su3(u, 1).re;
	hmc_float u_e21_im = reconstruct_su3(u, 1).im;
	hmc_float u_e22_re = reconstruct_su3(u, 2).re;
	hmc_float u_e22_im = reconstruct_su3(u, 2).im;
	
	tmp.e0.re = u.e00.re*in.e0.re + u.e01.re * in.e1.re + u.e02.re * in.e2.re
		  			- u.e00.im*in.e0.im - u.e01.im * in.e1.im - u.e02.im * in.e2.im;
	tmp.e0.im = u.e00.re*in.e0.im + u.e01.re * in.e1.im + u.e02.re * in.e2.im
		  			+ u.e00.im*in.e0.re + u.e01.im * in.e1.re + u.e02.im * in.e2.re;

	tmp.e1.re = u.e10.re*in.e0.re + u.e11.re * in.e1.re + u.e12.re * in.e2.re
		  			- u.e10.im*in.e0.im - u.e11.im * in.e1.im - u.e12.im * in.e2.im;
	tmp.e1.im = u.e10.re*in.e0.im + u.e11.re * in.e1.im + u.e12.re * in.e2.im
		  			+ u.e10.im*in.e0.re + u.e11.im * in.e1.re + u.e12.im * in.e2.re;

	tmp.e2.re = u_e20_re*in.e0.re + u_e21_re * in.e1.re + u_e22_re * in.e2.re
		  			- u_e20_im*in.e0.im - u_e21_im * in.e1.im - u_e22_im * in.e2.im;
	tmp.e2.im = u_e20_re*in.e0.im + u_e21_re * in.e1.im + u_e22_re * in.e2.im
		  			+ u_e20_im*in.e0.re + u_e21_im * in.e1.re + u_e22_im * in.e2.re;

	#else
	tmp.e0.re = u.e00.re*in.e0.re + u.e01.re * in.e1.re + u.e02.re * in.e2.re
		  			- u.e00.im*in.e0.im - u.e01.im * in.e1.im - u.e02.im * in.e2.im;
	tmp.e0.im = u.e00.re*in.e0.im + u.e01.re * in.e1.im + u.e02.re * in.e2.im
		  			+ u.e00.im*in.e0.re + u.e01.im * in.e1.re + u.e02.im * in.e2.re;

	tmp.e1.re = u.e10.re*in.e0.re + u.e11.re * in.e1.re + u.e12.re * in.e2.re
		  			- u.e10.im*in.e0.im - u.e11.im * in.e1.im - u.e12.im * in.e2.im;
	tmp.e1.im = u.e10.re*in.e0.im + u.e11.re * in.e1.im + u.e12.re * in.e2.im
		  			+ u.e10.im*in.e0.re + u.e11.im * in.e1.re + u.e12.im * in.e2.re;

	tmp.e2.re = u.e20.re*in.e0.re + u.e21.re * in.e1.re + u.e22.re * in.e2.re
		  			- u.e20.im*in.e0.im - u.e21.im * in.e1.im - u.e22.im * in.e2.im;
	tmp.e2.im = u.e20.re*in.e0.im + u.e21.re * in.e1.im + u.e22.re * in.e2.im
		  			+ u.e20.im*in.e0.re + u.e21.im * in.e1.re + u.e22.im * in.e2.re;

	#endif
	return tmp;
}

su3vec su3matrix_dagger_times_su3vec(Matrixsu3 u, su3vec in){
	su3vec tmp;
	#ifdef _RECONSTRUCT_TWELVE_
	hmc_float u_e20_re = reconstruct_su3(u, 0).re;
	hmc_float u_e20_im = reconstruct_su3(u, 0).im;
	hmc_float u_e21_re = reconstruct_su3(u, 1).re;
	hmc_float u_e21_im = reconstruct_su3(u, 1).im;
	hmc_float u_e22_re = reconstruct_su3(u, 2).re;
	hmc_float u_e22_im = reconstruct_su3(u, 2).im;
	
	tmp.e0.re = u.e00.re*in.e0.re + u.e10.re * in.e1.re + u_e20_re * in.e2.re
		  + u.e00.im*in.e0.im + u.e10.im * in.e1.im + u_e20_im * in.e2.im;
	tmp.e0.im = u.e00.re*in.e0.im + u.e10.re * in.e1.im + u_e20_re * in.e2.im
		  - u.e00.im*in.e0.re - u.e10.im * in.e1.re - u_e20_im * in.e2.re;

	tmp.e1.re = u.e01.re*in.e0.re + u.e11.re * in.e1.re + u_e21_re * in.e2.re
		  + u.e01.im*in.e0.im + u.e11.im * in.e1.im + u_e21_im * in.e2.im;
	tmp.e1.im = u.e01.re*in.e0.im + u.e11.re * in.e1.im + u_e21_re * in.e2.im
		  - u.e01.im*in.e0.re - u.e11.im * in.e1.re - u_e21_im * in.e2.re;

	tmp.e2.re = u.e02.re*in.e0.re + u.e12.re * in.e1.re + u_e22_re * in.e2.re
		  + u.e02.im*in.e0.im + u.e12.im * in.e1.im + u_e22_im * in.e2.im;
	tmp.e2.im = u.e02.re*in.e0.im + u.e12.re * in.e1.im + u_e22_re * in.e2.im
		  - u.e02.im*in.e0.re - u.e12.im * in.e1.re - u_e22_im * in.e2.re;

	#else
	tmp.e0.re = u.e00.re*in.e0.re + u.e10.re * in.e1.re + u.e20.re * in.e2.re
		  			+ u.e00.im*in.e0.im + u.e10.im * in.e1.im + u.e20.im * in.e2.im;
	tmp.e0.im = u.e00.re*in.e0.im + u.e10.re * in.e1.im + u.e20.re * in.e2.im
		 				- u.e00.im*in.e0.re - u.e10.im * in.e1.re - u.e20.im * in.e2.re;

	tmp.e1.re = u.e01.re*in.e0.re + u.e11.re * in.e1.re + u.e21.re * in.e2.re
		  			+ u.e01.im*in.e0.im + u.e11.im * in.e1.im + u.e21.im * in.e2.im;
	tmp.e1.im = u.e01.re*in.e0.im + u.e11.re * in.e1.im + u.e21.re * in.e2.im
		  			- u.e01.im*in.e0.re - u.e11.im * in.e1.re - u.e21.im * in.e2.re;

	tmp.e2.re = u.e02.re*in.e0.re + u.e12.re * in.e1.re + u.e22.re * in.e2.re
		  			+ u.e02.im*in.e0.im + u.e12.im * in.e1.im + u.e22.im * in.e2.im;
	tmp.e2.im = u.e02.re*in.e0.im + u.e12.re * in.e1.im + u.e22.re * in.e2.im
		  			- u.e02.im*in.e0.re - u.e12.im * in.e1.re - u.e22.im * in.e2.re;

	#endif


	return tmp;
}

su3vec su3vec_times_minusone(su3vec in){
	su3vec tmp;
	tmp.e0.re = -(in).e0.re;
	tmp.e0.im = -(in).e0.im;
	tmp.e1.re = -(in).e1.re;
	tmp.e1.im = -(in).e1.im;
	tmp.e2.re = -(in).e2.re;
	tmp.e2.im = -(in).e2.im;
	return tmp;
}

su3vec su3vec_acc(su3vec in1, su3vec in2){
	su3vec tmp;     
	tmp.e0.re = in1.e0.re + in2.e0.re;
	tmp.e0.im = in1.e0.im + in2.e0.im;
	tmp.e1.re = in1.e1.re + in2.e1.re;
	tmp.e1.im = in1.e1.im + in2.e1.im;
	tmp.e2.re = in1.e2.re + in2.e2.re;
	tmp.e2.im = in1.e2.im + in2.e2.im;
	return tmp;
}

su3vec su3vec_acc_i(su3vec in1, su3vec in2){
	su3vec tmp;     
	tmp.e0.re = in1.e0.re - in2.e0.im;
	tmp.e0.im = in1.e0.im + in2.e0.re;
	tmp.e1.re = in1.e1.re - in2.e1.im;
	tmp.e1.im = in1.e1.im + in2.e1.re;
	tmp.e2.re = in1.e2.re - in2.e2.im;
	tmp.e2.im = in1.e2.im + in2.e2.re;
	return tmp;
}

su3vec su3vec_dim(su3vec in1, su3vec in2){
	su3vec tmp;     
	tmp.e0.re = in1.e0.re - in2.e0.re;
	tmp.e0.im = in1.e0.im - in2.e0.im;
	tmp.e1.re = in1.e1.re - in2.e1.re;
	tmp.e1.im = in1.e1.im - in2.e1.im;
	tmp.e2.re = in1.e2.re - in2.e2.re;
	tmp.e2.im = in1.e2.im - in2.e2.im;
	return tmp;
}

su3vec su3vec_dim_i(su3vec in1, su3vec in2){
	su3vec tmp;     
	tmp.e0.re = in1.e0.re + in2.e0.im;
	tmp.e0.im = in1.e0.im - in2.e0.re;
	tmp.e1.re = in1.e1.re + in2.e1.im;
	tmp.e1.im = in1.e1.im - in2.e1.re;
	tmp.e2.re = in1.e2.re + in2.e2.im;
	tmp.e2.im = in1.e2.im - in2.e2.re;
	return tmp;
}

su3vec su3vec_acc_acc(su3vec in1, su3vec in2, su3vec in3){
	su3vec tmp;     
	tmp.e0.re = in1.e0.re + in2.e0.re + in3.e0.re;
	tmp.e0.im = in1.e0.im + in2.e0.im + in3.e0.im;
	tmp.e1.re = in1.e1.re + in2.e1.re + in3.e1.re;
	tmp.e1.im = in1.e1.im + in2.e1.im + in3.e1.im;
	tmp.e2.re = in1.e2.re + in2.e2.re + in3.e2.re;
	tmp.e2.im = in1.e2.im + in2.e2.im + in3.e2.im;
	return tmp;
}


//calculates the Dirac-Trace of the matrix resulting from multiplying U*V^dagger =  u*v^dagger + w*x^dagger, where u, v, w, x are SU(3)-vectors (using spinprojection). The result is a 3x3-matrix
Matrix3x3 tr_v_times_u_dagger(su3vec u, su3vec v, su3vec w, su3vec x){
	Matrix3x3 tmp;
	tmp.e00.re=(u).e0.re*(v).e0.re+(u).e0.im*(v).e0.im + (w).e0.re*(x).e0.re+(w).e0.im*(x).e0.im; 
	tmp.e00.im=(u).e0.im*(v).e0.re-(u).e0.re*(v).e0.im + (w).e0.im*(x).e0.re-(w).e0.re*(x).e0.im;
	tmp.e01.re=(u).e0.re*(v).e1.re+(u).e0.im*(v).e1.im + (w).e0.re*(x).e1.re+(w).e0.im*(x).e1.im; 
	tmp.e01.im=(u).e0.im*(v).e1.re-(u).e0.re*(v).e1.im + (w).e0.im*(x).e1.re-(w).e0.re*(x).e1.im;
	tmp.e02.re=(u).e0.re*(v).e2.re+(u).e0.im*(v).e2.im + (w).e0.re*(x).e2.re+(w).e0.im*(x).e2.im; 
	tmp.e02.im=(u).e0.im*(v).e2.re-(u).e0.re*(v).e2.im + (w).e0.im*(x).e2.re-(w).e0.re*(x).e2.im; 
	tmp.e10.re=(u).e1.re*(v).e0.re+(u).e1.im*(v).e0.im + (w).e1.re*(x).e0.re+(w).e1.im*(x).e0.im; 
	tmp.e10.im=(u).e1.im*(v).e0.re-(u).e1.re*(v).e0.im + (w).e1.im*(x).e0.re-(w).e1.re*(x).e0.im; 
	tmp.e11.re=(u).e1.re*(v).e1.re+(u).e1.im*(v).e1.im + (w).e1.re*(x).e1.re+(w).e1.im*(x).e1.im; 
	tmp.e11.im=(u).e1.im*(v).e1.re-(u).e1.re*(v).e1.im + (w).e1.im*(x).e1.re-(w).e1.re*(x).e1.im; 
	tmp.e12.re=(u).e1.re*(v).e2.re+(u).e1.im*(v).e2.im + (w).e1.re*(x).e2.re+(w).e1.im*(x).e2.im; 
	tmp.e12.im=(u).e1.im*(v).e2.re-(u).e1.re*(v).e2.im + (w).e1.im*(x).e2.re-(w).e1.re*(x).e2.im; 
	tmp.e20.re=(u).e2.re*(v).e0.re+(u).e2.im*(v).e0.im + (w).e2.re*(x).e0.re+(w).e2.im*(x).e0.im; 
	tmp.e20.im=(u).e2.im*(v).e0.re-(u).e2.re*(v).e0.im + (w).e2.im*(x).e0.re-(w).e2.re*(x).e0.im; 
	tmp.e21.re=(u).e2.re*(v).e1.re+(u).e2.im*(v).e1.im + (w).e2.re*(x).e1.re+(w).e2.im*(x).e1.im; 
	tmp.e21.im=(u).e2.im*(v).e1.re-(u).e2.re*(v).e1.im + (w).e2.im*(x).e1.re-(w).e2.re*(x).e1.im; 
	tmp.e22.re=(u).e2.re*(v).e2.re+(u).e2.im*(v).e2.im + (w).e2.re*(x).e2.re+(w).e2.im*(x).e2.im; 
	tmp.e22.im=(u).e2.im*(v).e2.re-(u).e2.re*(v).e2.im + (w).e2.im*(x).e2.re-(w).e2.re*(x).e2.im; 
	return tmp;
}
/** @file
 * Device code for operations on spinors
 */

//opencl_operations_spinor

#ifdef _FERMIONS_
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
	hmc_float res=0;
	res += su3vec_squarenorm(in.e0);
	res += su3vec_squarenorm(in.e1);
	res += su3vec_squarenorm(in.e2);
	res += su3vec_squarenorm(in.e3);
	return res;
}

hmc_complex spinor_scalarproduct(spinor in1, spinor in2){
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



//CP: I dont think this function is needed anywhere?!?!
// void spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr)
// {
// 	for(int j=0; j<SPINORSIZE; j++) {
// 		inout[j].re += incr[j].re;
// 		inout[j].im += incr[j].im;
// 	}
// 	return;
// }

//CP: I think this can be done by using su3vec_times_complex...
// void spinor_apply_bc(hmc_spinor * in, hmc_float theta)
// {
// 	for(int n = 0; n<SPINORSIZE; n++) {
// 		hmc_float tmp1 = in[n].re;
// 		hmc_float tmp2 = in[n].im;
// 		in[n].re = cos(theta)*tmp1 - sin(theta)*tmp2;
// 		in[n].im = sin(theta)*tmp1 + cos(theta)*tmp2;
// 	}
// 	return;
// }

//CP: I think this is not needed anymore. A more flexibel function will be introduced that can be used for the inverse sitediagonal matrix as well
/*
//spinout =  (1 + 2*i*gamma_5*kappa*mu)spin_in
spinor M_diag_local(spinor y, hmc_float mubar)
{
	spinor out_tmp;

	//CP: how is this called now??
	#ifdef _TWISTEDMASS_
	hmc_complex twistfactor = {1., mubar};
	hmc_complex twistfactor_minus = {1., -1.*mubar};
	//Diagonalpart:
	//	(1+i*mubar*gamma_5)psi = (1, mubar)psi.0,1 (1,-mubar)psi.2,3
        out_tmp.e0 = su3vec_times_complex(y.e0, twistfactor);
	out_tmp.e1 = su3vec_times_complex(y.e1, twistfactor);
	out_tmp.e2 = su3vec_times_complex(y.e2, twistfactor_minus);
	out_tmp.e3 = su3vec_times_complex(y.e3, twistfactor_minus);
	#else
	//Pure Wilson:
	//Diagonalpart:	(1)psi -> spinor not changed!!
	#endif
	return out_tmp;
}
*/

//CP: this is all deprecated!!!
/*
//spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
void dslash_0(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma0(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma0(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}

// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
void dslash_1(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma1(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma1(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}

// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
void dslash_2(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma2(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma2(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}

// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
void dslash_3(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_ocl_su3matrix* u, hmc_ocl_su3matrix* udagger)
{

	spinprojectproduct_gamma3(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma3(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);

	return;
}
*/
#endif //_FERMIONS_
/** @file
 * Device code for operations on the spinor field
 */


void print_su3vec(su3vec in){
     printf("(%f,%f)\t(%f,%f)\t(%f,%f)\t", in.e0.re, in.e0.im, in.e1.re, in.e1.im, in.e2.re, in.e2.im);
}

void print_spinor(spinor in){
     print_su3vec(in.e0);
     print_su3vec(in.e1);
     print_su3vec(in.e2);
     print_su3vec(in.e3);
     printf("\n");
}


spinor get_spinor_from_field(__global spinorfield* in, int n, int t)
{
	int pos = get_global_pos(n,t);
	spinor out;
	out = in[pos];
	return out;
}

void put_spinor_to_field(spinor in, __global spinorfield* out, int n, int t)
{
	int pos = get_global_pos(n,t);
	out[pos] = in;
}
/**
 @file fermionmatrix-functions
*/

//local twisted-mass Diagonalmatrix:
//	(1+i*mubar*gamma_5)psi = (1, mubar)psi.0,1 (1,-mubar)psi.2,3
spinor inline M_diag_tm_local(spinor const in, hmc_complex const factor1, hmc_complex const factor2){
	spinor tmp;
	tmp.e0 = su3vec_times_complex(in.e0, factor1);
	tmp.e1 = su3vec_times_complex(in.e1, factor1);
	tmp.e2 = su3vec_times_complex(in.e2, factor2);
	tmp.e3 = su3vec_times_complex(in.e3, factor2);
	return tmp;	
}

/** @todo this can be optimized... */
//local gamma5:
//	(gamma_5)psi = (1)psi.0,1 (-1)psi.2,3
spinor inline gamma5_local(spinor const in){
	spinor tmp;
	tmp.e0 = in.e0;
	tmp.e1 = in.e1;
	tmp.e2 = su3vec_times_real(in.e2, -1.);
	tmp.e3 = su3vec_times_real(in.e3, -1.);
	return tmp;	
}

//"local" dslash working on a particular link (n,t) in a specific direction
//NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the correct dslash!!!
//spinor dslash_local_0(__global const spinorfield * const restrict in,__global const ocl_s_gaugefield * const restrict field, const int n, const int t){
spinor dslash_local_0(__global spinorfield * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();


	//go through the different directions
	///////////////////////////////////
	// mu = 0
	///////////////////////////////////
	dir = 0;
	///////////////////////////////////
	//mu = +0
	nn = (t+1)%NTIME;
	plus = get_spinor_from_field(in, n, nn);
	U = field[get_global_link_pos(dir, n, t)];
	//if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
	bc_tmp.re = KAPPA_TEMPORAL_RE;
	bc_tmp.im = KAPPA_TEMPORAL_IM;
	///////////////////////////////////
	// Calculate psi/phi = (1 - gamma_0) plus/y
	// with 1 - gamma_0:
	// | 1  0  1  0 |        | psi.e0 + psi.e2 |
	// | 0  1  0  1 |  psi = | psi.e1 + psi.e3 |
	// | 1  0  1  0 |        | psi.e1 + psi.e3 |
	// | 0  1  0  1 |        | psi.e0 + psi.e2 |
	///////////////////////////////////
	// psi = 0. component of (1-gamma_0)y
	psi = su3vec_acc(plus.e0, plus.e2);
	// phi = U*psi
	phi =  su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);
	// psi = 1. component of (1-gamma_0)y
	psi = su3vec_acc(plus.e1, plus.e3);
	// phi = U*psi
	phi =  su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);

	/////////////////////////////////////
	//mu = -0
	nn = (t-1+NTIME)%NTIME;
	plus = get_spinor_from_field(in, n, nn);
	U = field[get_global_link_pos(dir, n, nn)];
	//if chemical potential is activated, U has to be multiplied by appropiate factor
	//this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
	//	(exp(iq)U)^+ = exp(-iq)U^+
	//as it should be
	//in the real case, one has to take exp(q) -> exp(-q)
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, MEXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_TEMPORAL_RE;
	bc_tmp.im = MKAPPA_TEMPORAL_IM;
	///////////////////////////////////
	// Calculate psi/phi = (1 + gamma_0) y
	// with 1 + gamma_0:
	// | 1  0 -1  0 |       | psi.e0 - psi.e2 |
	// | 0  1  0 -1 | psi = | psi.e1 - psi.e3 |
	// |-1  0  1  0 |       | psi.e1 - psi.e2 |
	// | 0 -1  0  1 |       | psi.e0 - psi.e3 |
	///////////////////////////////////
	// psi = 0. component of (1+gamma_0)y
	psi = su3vec_dim(plus.e0, plus.e2);
	// phi = U*psi
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);
	// psi = 1. component of (1+gamma_0)y
	psi = su3vec_dim(plus.e1, plus.e3);
	// phi = U*psi
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);		

	return out_tmp;
}

spinor dslash_local_1(__global spinorfield * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();

	//CP: all actions correspond to the mu = 0 ones
	///////////////////////////////////
	// mu = 1
	///////////////////////////////////
	dir = 1;
	
	///////////////////////////////////
	// mu = +1
	nn = get_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
	U = field[get_global_link_pos(dir, n, t)];
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = KAPPA_SPATIAL_IM;
	/////////////////////////////////
	//Calculate (1 - gamma_1) y
	//with 1 - gamma_1:
	//| 1  0  0  i |       |       psi.e0 + i*psi.e3  |
	//| 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
	//| 0  i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
	//| i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
	/////////////////////////////////
	psi = su3vec_acc_i(plus.e0, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);
	
	psi = su3vec_acc_i(plus.e1, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);		

	///////////////////////////////////
	//mu = -1
	nn = get_lower_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
	U = field[get_global_link_pos(dir, nn, t)];
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = MKAPPA_SPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 + gamma_1) y
	// with 1 + gamma_1:
	// | 1  0  0 -i |       |       psi.e0 - i*psi.e3  |
	// | 0  1 -i  0 | psi = |       psi.e1 - i*psi.e2  |
	// | 0  i  1  0 |       |(-i)*( psi.e1 - i*psi.e2) |
	// | i  0  0  1 |       |(-i)*( psi.e0 - i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim_i(plus.e0, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);
	
	psi = su3vec_dim_i(plus.e1, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);
	
	
	return out_tmp;
}

spinor dslash_local_2(__global spinorfield * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();;
	
	///////////////////////////////////
	// mu = 2
	///////////////////////////////////
	dir = 2;
	
	///////////////////////////////////
	// mu = +2
	nn = get_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
	U = field[get_global_link_pos(dir, n, t)];
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = KAPPA_SPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 - gamma_2) y
	// with 1 - gamma_2:
	// | 1  0  0  1 |       |       psi.e0 + psi.e3  |
	// | 0  1 -1  0 | psi = |       psi.e1 - psi.e2  |
	// | 0 -1  1  0 |       |(-1)*( psi.e1 + psi.e2) |
	// | 1  0  0  1 |       |     ( psi.e0 + psi.e3) |
	///////////////////////////////////
	psi = su3vec_acc(plus.e0, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);
	
	psi = su3vec_dim(plus.e1, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);	
	
	///////////////////////////////////
	//mu = -2
	nn = get_lower_neighbor(n,dir);
	plus = get_spinor_from_field(in,  nn, t);
	U = field[get_global_link_pos(dir, nn, t)];
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = MKAPPA_SPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 + gamma_2) y
	// with 1 + gamma_2:
	// | 1  0  0 -1 |       |       psi.e0 - psi.e3  |
	// | 0  1  1  0 | psi = |       psi.e1 + psi.e2  |
	// | 0  1  1  0 |       |     ( psi.e1 + psi.e2) |
	// |-1  0  0  1 |       |(-1)*( psi.e0 - psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim(plus.e0, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);
	
	psi = su3vec_acc(plus.e1, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);
	
	return out_tmp;
}

spinor dslash_local_3(__global spinorfield * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();
	
	///////////////////////////////////
	// mu = 3
	///////////////////////////////////
	dir = 3;
	
	///////////////////////////////////
	// mu = +3
	nn = get_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
	U = field[get_global_link_pos(dir, n, t)];
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = KAPPA_SPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 - gamma_3) y
	// with 1 - gamma_3:
	// | 1  0  i  0 |        |       psi.e0 + i*psi.e2  |
	// | 0  1  0 -i |  psi = |       psi.e1 - i*psi.e3  |
	// |-i  0  1  0 |        |   i *(psi.e0 + i*psi.e2) |
	// | 0  i  0  1 |        | (-i)*(psi.e1 - i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_acc_i(plus.e0, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);
	
	psi = su3vec_dim_i(plus.e1, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);	

	///////////////////////////////////
	//mu = -3
	nn = get_lower_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
	U = field[get_global_link_pos(dir, nn, t)];
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = MKAPPA_SPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 + gamma_3) y
	// with 1 + gamma_3:
	// | 1  0 -i  0 |       |       psi.e0 - i*psi.e2  |
	// | 0  1  0  i | psi = |       psi.e1 + i*psi.e3  |
	// | i  0  1  0 |       | (-i)*(psi.e0 - i*psi.e2) |
	// | 0 -i  0  1 |       |   i *(psi.e1 + i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim_i(plus.e0, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);
	
	psi = su3vec_acc_i(plus.e1, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);
	
	return out_tmp;
}

__kernel void M_tm_plus(__global spinorfield * in,  __global ocl_s_gaugefield * field, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp, out_tmp2;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);

		//get input spinor
		plus = get_spinor_from_field(in, n, t);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
		out_tmp2 = set_spinor_zero();

		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local_0(in, field, n, t);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_1(in, field, n, t);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_2(in, field, n, t);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_3(in, field, n, t);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		put_spinor_to_field(out_tmp, out, n, t);
	}
}
