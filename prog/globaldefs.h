/** @file
 * Global parameters, specified as macros
 *
 * @todo Some, or even most of these, would be better represented as constants
 */

#ifndef _GLOBALSH_
#define _GLOBALSH_

/** Number of colors */
#define NC 3
#define NSPIN 4

/** Number of dimensions of the lattice */
#define NDIM 4

#ifndef _INKERNEL_

/** Spatial volume of the lattice */
#define VOLSPACE NSPACE*NSPACE*NSPACE
/** 4-Dimensional Volume of the lattice */
#define VOL4D VOLSPACE*NTIME

#define SPINORSIZE NSPIN*NC
#define HALFSPINORSIZE NSPIN/2*NC
//#define SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE
#define SPINORFIELDSIZE NTIME*VOLSPACE
#define EOPREC_SPINORFIELDSIZE NTIME*VOLSPACE/2
//CP: This is not needed anywhere!!
//#define GAUGEMOMENTASIZE NDIM*VOL4D*(NC*NC-1)
#define GAUGEFIELDSIZE NC*NC*NDIM*VOL4D
#define SU3ALGEBRASIZE NC*NC-1

#define GAUGEMOMENTASIZE2 NDIM*VOL4D

//startconditions:
#define START_FROM_SOURCE 2
#define COLD_START 0
#define HOT_START 1

#endif //_INKERNEL_

#define WILSON 0
#define CLOVER 1
#define TWISTEDMASS 2

//EVEN ODD
#define EVEN 0
#define ODD 1

/**
 * PI
 * @todo Rather use PI from stdlib
 */
#define PI 	3.14159265358979

#define su2_entries 4

/** Number of threads to use for OpenCL kernels */
#ifdef _USEGPU_
#define NUMTHREADS 128
#else
#define NUMTHREADS 1
#endif

#ifndef _INKERNEL_
// Definition of numeric constants for the symmetric structure constants d_ijk of su(3)
/** 1/2 */
#define F_1_2   (static_cast<hmc_float>(0.5))
/** 1/(2*sqrt(3)) */
#define F_1_2S3 (static_cast<hmc_float>(0.288675134594813))
/** 1/sqrt(3) */
#define F_1_S3  (static_cast<hmc_float>(0.577350269189626))

#endif //_INKERNEL_

#endif


