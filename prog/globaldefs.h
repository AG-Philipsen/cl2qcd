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

//EVEN ODD
#define EVEN 0
#define ODD 1
#define OE 0
#define EO 1


#define ARG_DEF -1

/**
 * PI
 * @todo Rather use PI from stdlib
 */
#define PI  3.14159265358979

#define su2_entries 4

// Definition of numeric constants for the symmetric structure constants d_ijk of su(3) suited for OpenCL
#ifdef _INKERNEL_
/** 1/2 */
#define F_1_2  0.5
/** 1/(2*sqrt(3)) */
#define F_1_2S3 0.288675134594813
/** 1/sqrt(3) */
#define F_1_S3  0.577350269189626
#else
/** 1/2 */
#define F_1_2   (static_cast<hmc_float>(0.5))
/** 1/(2*sqrt(3)) */
#define F_1_2S3 (static_cast<hmc_float>(0.288675134594813))
/** 1/sqrt(3) */
#define F_1_S3  (static_cast<hmc_float>(0.577350269189626))
#endif //_INKERNEL_

#endif


