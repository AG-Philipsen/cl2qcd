#include "globaldefs.h" //NDIM, NSPIN, NC
#include "types.h"

//for hmc_ocl_su3matrix
#ifdef _RECONSTRUCT_TWELVE_
#define SU3SIZE NC*(NC-1)
#else
#define SU3SIZE NC*NC
#endif

