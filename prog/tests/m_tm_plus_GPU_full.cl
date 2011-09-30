/** @file
 * All definitions, types and functions required for the fermionmatrix M_tm_plus
 * No more files are included, everything that is not used has been deleted.
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

//globaldefs.h

/** Number of colors */
#define NC 3
#define NSPIN 4

/** Number of dimensions of the lattice */
#define NDIM 4

#ifndef _INKERNEL_

#define SPINORSIZE NSPIN*NC
#define HALFSPINORSIZE NSPIN/2*NC
#define SU3ALGEBRASIZE NC*NC-1

//startconditions:
#define START_FROM_SOURCE 2
#define COLD_START 0
#define HOT_START 1

#endif //_INKERNEL_

#define WILSON 0
#define CLOVER 1
#define TWISTEDMASS 2

#define LEAPFROG 0
#define TWOMN 1

//EVEN ODD
#define EVEN 0
#define ODD 1

/**
 * PI
 * @todo Rather use PI from stdlib
 */
#define PI  3.14159265358979

#define su2_entries 4

//types.h

#ifndef _INKERNEL_
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#endif

#ifdef _INKERNEL_
#define CONST __constant
#else /* _INKERNEL_ */
#define CONST const
#endif /* _INKERNEL_ */

/** The floating precision type used by hmc, can be 32 or 64 bit. */
#ifdef _USEDOUBLEPREC_
typedef double hmc_float;
#else
typedef float hmc_float;
#endif

/** An OpenCL-compatible constant 1.
 * @todo this is rediculous, 1.0(f) is way more readable, let the
 *       compiler do the optimizations.
 */
#ifdef _INKERNEL_
__constant hmc_float hmc_one_f = 1.0f;
#else
hmc_float const hmc_one_f = static_cast<hmc_float>(1);
#endif

/** Complex number type, precision is the same as for hmc_float */
#ifdef _INKERNEL_
typedef struct {
	hmc_float re;
	hmc_float im;
} hmc_complex;
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e10;
	hmc_complex e11;
} Matrixsu2;
typedef struct {
	hmc_float e00;
	hmc_float e01;
	hmc_float e10;
	hmc_float e11;
} Matrixsu2_pauli;
#else
struct hmc_complex {
	hmc_float re;
	hmc_float im;
};
struct Matrixsu2 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e10;
	hmc_complex e11;
} ;
struct Matrixsu2_pauli {
	hmc_float e00;
	hmc_float e01;
	hmc_float e10;
	hmc_float e11;
} ;
#endif

#ifdef _INKERNEL_
__constant hmc_complex hmc_complex_one = {1., 0.};
__constant hmc_complex hmc_complex_zero = {0., 0.};
__constant hmc_complex hmc_complex_minusone = { -1., 0.};
__constant hmc_complex hmc_complex_i = {0., 1.};
#else
/** A complex 1 */
hmc_complex const hmc_complex_one = {1., 0.};
/** A complex -1 */
hmc_complex const hmc_complex_minusone = { -1., 0.};
/** A complex 0 */
hmc_complex const hmc_complex_zero = {0., 0.};
/** A complex i */
hmc_complex const hmc_complex_i = {0., 1.};
#endif

//also CPU
#ifndef _INKERNEL_

//define a gauge field: gauge[su3][mu][coord3d][coord_time]
#ifdef _RECONSTRUCT_TWELVE_
typedef hmc_complex hmc_su3matrix [NC*(NC-1)];

struct Matrixsu3 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
};

typedef hmc_complex hmc_staplematrix [NC*NC];
typedef hmc_complex hmc_3x3matrix[3][3];
#else //_RECONSTRUCT_TWELVE_
/** A generic SU3 matrix */
typedef hmc_complex hmc_su3matrix [NC][NC];
/** A matrix representing a staple */
typedef hmc_su3matrix hmc_staplematrix;
/** A generic 3x3 matrix */
typedef hmc_complex hmc_3x3matrix[3][3];

struct Matrixsu3 {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
};
#endif //_RECONSTRUCT_TWELVE_

typedef Matrixsu3 ocl_s_gaugefield;

#endif // ifndef _INKERNEL_

typedef hmc_complex hmc_ocl_su3matrix;
typedef hmc_complex hmc_ocl_3x3matrix;
typedef hmc_complex hmc_ocl_staplematrix;
typedef hmc_float hmc_ocl_gaugefield;

#ifdef _INKERNEL_

typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
} Matrix3x3;

#ifdef _RECONSTRUCT_TWELVE_
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
} Matrixsu3;
#else
typedef struct {
	hmc_complex e00;
	hmc_complex e01;
	hmc_complex e02;
	hmc_complex e10;
	hmc_complex e11;
	hmc_complex e12;
	hmc_complex e20;
	hmc_complex e21;
	hmc_complex e22;
} Matrixsu3;

#endif //ifdef _REC12_

#endif  //ifdef _INKERNEL_

typedef Matrixsu3 ocl_s_gaugefield;

#ifndef _INKERNEL_
/** Storage type for state of the random number generator */
typedef cl_ulong4 hmc_ocl_ran;
#endif /* _INKERNEL_ */

//CP: this an algebraelement
typedef struct {
	hmc_float e0;
	hmc_float e1;
	hmc_float e2;
	hmc_float e3;
	hmc_float e4;
	hmc_float e5;
	hmc_float e6;
	hmc_float e7;
} ae;

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
 * Common fermion types used by HMC, both on host and OpenCL device.
 */

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

// su3vec su3vec_times_minusone(su3vec in){
// 	su3vec tmp;
// 	tmp.e0.re = -(in).e0.re;
// 	tmp.e0.im = -(in).e0.im;
// 	tmp.e1.re = -(in).e1.re;
// 	tmp.e1.im = -(in).e1.im;
// 	tmp.e2.re = -(in).e2.re;
// 	tmp.e2.im = -(in).e2.im;
// 	return tmp;
// }

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

// su3vec su3vec_acc_acc(su3vec in1, su3vec in2, su3vec in3){
// 	su3vec tmp;     
// 	tmp.e0.re = in1.e0.re + in2.e0.re + in3.e0.re;
// 	tmp.e0.im = in1.e0.im + in2.e0.im + in3.e0.im;
// 	tmp.e1.re = in1.e1.re + in2.e1.re + in3.e1.re;
// 	tmp.e1.im = in1.e1.im + in2.e1.im + in3.e1.im;
// 	tmp.e2.re = in1.e2.re + in2.e2.re + in3.e2.re;
// 	tmp.e2.im = in1.e2.im + in2.e2.im + in3.e2.im;
// 	return tmp;
// }


// //calculates the Dirac-Trace of the matrix resulting from multiplying U*V^dagger =  u*v^dagger + w*x^dagger, where u, v, w, x are SU(3)-vectors (using spinprojection). The result is a 3x3-matrix
// Matrix3x3 tr_v_times_u_dagger(su3vec u, su3vec v, su3vec w, su3vec x){
// 	Matrix3x3 tmp;
// 	tmp.e00.re=(u).e0.re*(v).e0.re+(u).e0.im*(v).e0.im + (w).e0.re*(x).e0.re+(w).e0.im*(x).e0.im; 
// 	tmp.e00.im=(u).e0.im*(v).e0.re-(u).e0.re*(v).e0.im + (w).e0.im*(x).e0.re-(w).e0.re*(x).e0.im;
// 	tmp.e01.re=(u).e0.re*(v).e1.re+(u).e0.im*(v).e1.im + (w).e0.re*(x).e1.re+(w).e0.im*(x).e1.im; 
// 	tmp.e01.im=(u).e0.im*(v).e1.re-(u).e0.re*(v).e1.im + (w).e0.im*(x).e1.re-(w).e0.re*(x).e1.im;
// 	tmp.e02.re=(u).e0.re*(v).e2.re+(u).e0.im*(v).e2.im + (w).e0.re*(x).e2.re+(w).e0.im*(x).e2.im; 
// 	tmp.e02.im=(u).e0.im*(v).e2.re-(u).e0.re*(v).e2.im + (w).e0.im*(x).e2.re-(w).e0.re*(x).e2.im; 
// 	tmp.e10.re=(u).e1.re*(v).e0.re+(u).e1.im*(v).e0.im + (w).e1.re*(x).e0.re+(w).e1.im*(x).e0.im; 
// 	tmp.e10.im=(u).e1.im*(v).e0.re-(u).e1.re*(v).e0.im + (w).e1.im*(x).e0.re-(w).e1.re*(x).e0.im; 
// 	tmp.e11.re=(u).e1.re*(v).e1.re+(u).e1.im*(v).e1.im + (w).e1.re*(x).e1.re+(w).e1.im*(x).e1.im; 
// 	tmp.e11.im=(u).e1.im*(v).e1.re-(u).e1.re*(v).e1.im + (w).e1.im*(x).e1.re-(w).e1.re*(x).e1.im; 
// 	tmp.e12.re=(u).e1.re*(v).e2.re+(u).e1.im*(v).e2.im + (w).e1.re*(x).e2.re+(w).e1.im*(x).e2.im; 
// 	tmp.e12.im=(u).e1.im*(v).e2.re-(u).e1.re*(v).e2.im + (w).e1.im*(x).e2.re-(w).e1.re*(x).e2.im; 
// 	tmp.e20.re=(u).e2.re*(v).e0.re+(u).e2.im*(v).e0.im + (w).e2.re*(x).e0.re+(w).e2.im*(x).e0.im; 
// 	tmp.e20.im=(u).e2.im*(v).e0.re-(u).e2.re*(v).e0.im + (w).e2.im*(x).e0.re-(w).e2.re*(x).e0.im; 
// 	tmp.e21.re=(u).e2.re*(v).e1.re+(u).e2.im*(v).e1.im + (w).e2.re*(x).e1.re+(w).e2.im*(x).e1.im; 
// 	tmp.e21.im=(u).e2.im*(v).e1.re-(u).e2.re*(v).e1.im + (w).e2.im*(x).e1.re-(w).e2.re*(x).e1.im; 
// 	tmp.e22.re=(u).e2.re*(v).e2.re+(u).e2.im*(v).e2.im + (w).e2.re*(x).e2.re+(w).e2.im*(x).e2.im; 
// 	tmp.e22.im=(u).e2.im*(v).e2.re-(u).e2.re*(v).e2.im + (w).e2.im*(x).e2.re-(w).e2.re*(x).e2.im; 
// 	return tmp;
// }
/** @file
 * Device code for operations on spinors
 */

//opencl_operations_spinor

spinor set_spinor_zero()
{
	spinor tmp;
	tmp.e0 = set_su3vec_zero();
	tmp.e1 = set_su3vec_zero();
	tmp.e2 = set_su3vec_zero();
	tmp.e3 = set_su3vec_zero();
	return tmp;
}

// spinor set_spinor_cold()
// {
// 	spinor tmp;
// 	tmp.e0 = set_su3vec_cold();
// 	tmp.e1 = set_su3vec_cold();
// 	tmp.e2 = set_su3vec_cold();
// 	tmp.e3 = set_su3vec_cold();
// 	return tmp;
// }
 
// hmc_float spinor_squarenorm(spinor in)
// {
// 	hmc_float res=0;
// 	res += su3vec_squarenorm(in.e0);
// 	res += su3vec_squarenorm(in.e1);
// 	res += su3vec_squarenorm(in.e2);
// 	res += su3vec_squarenorm(in.e3);
// 	return res;
// }

// hmc_complex spinor_scalarproduct(spinor in1, spinor in2){
// 	hmc_complex res = hmc_complex_zero;
// 	hmc_complex tmp;
// 	tmp = su3vec_scalarproduct(in1.e0, in2.e0);
// 	res.re += tmp.re;
// 	res.im += tmp.im;
// 	tmp = su3vec_scalarproduct(in1.e1, in2.e1);
// 	res.re += tmp.re;
// 	res.im += tmp.im;
// 	tmp = su3vec_scalarproduct(in1.e2, in2.e2);
// 	res.re += tmp.re;
// 	res.im += tmp.im;
// 	tmp = su3vec_scalarproduct(in1.e3, in2.e3);
// 	res.re += tmp.re;
// 	res.im += tmp.im;
// 	return res;
// }

// spinor real_multiply_spinor(spinor in, hmc_float factor)
// {
// 	spinor tmp;
// 	tmp.e0 = su3vec_times_real(in.e0, factor);
// 	tmp.e1 = su3vec_times_real(in.e1, factor);
// 	tmp.e2 = su3vec_times_real(in.e2, factor);
// 	tmp.e3 = su3vec_times_real(in.e3, factor);
// 	return tmp;
// }

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

// spinor spinor_acc_acc(spinor in1, spinor in2, spinor in3)
// {
// 	spinor tmp;
// 	tmp.e0 = su3vec_acc_acc(in1.e0, in2.e0, in3.e0);
// 	tmp.e1 = su3vec_acc_acc(in1.e1, in2.e1, in3.e1);
// 	tmp.e2 = su3vec_acc_acc(in1.e2, in2.e2, in3.e2);
// 	tmp.e3 = su3vec_acc_acc(in1.e3, in2.e3, in3.e3);
// 	return tmp;
// }

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
