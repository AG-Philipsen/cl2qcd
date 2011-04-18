
#ifndef _GLOBALSH_
#define _GLOBALSH_

#define NC 3
#define NSPIN 4
#define NDIM 4

#ifndef _INKERNEL_


#define NTIME 4
#define NSPACE 8

#define VOLSPACE NSPACE*NSPACE*NSPACE
#define VOL4D VOLSPACE*NTIME

#define SPINORSIZE NSPIN*NC
#define HALFSPINORSIZE NSPIN/2*NC
#define SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE
#define EOPREC_SPINORFIELDSIZE SPINORSIZE*NTIME*VOLSPACE/2
#define GAUGEMOMENTASIZE NDIM*VOL4D*(NC*NC-1)
#define GAUGEFIELDSIZE NC*NC*NDIM*VOL4D

//startconditions:
#define START_FROM_SOURCE 2
#define COLD_START 0
#define HOT_START 1

#endif //_INKERNEL_

//EVEN ODD
#define EVEN 0
#define ODD 1

#define TRUE 1
#define FALSE 0

#define PI 	3.14159265358979

#define su2_entries 4

#ifdef _USEGPU_
#define NUMTHREADS 128
#else
#define NUMTHREADS 1
#endif

// Definition of numeric constants for the symmetric structure constants d_ijk of su(3)
#define F_1_2   (static_cast<hmc_float>(0.5))					// 1/2
#define F_1_2S3 (static_cast<hmc_float>(0.288675134594813))		// 1/(2*sqrt(3))
#define F_1_S3  (static_cast<hmc_float>(0.577350269189626))		// 1/sqrt(3)


#endif


