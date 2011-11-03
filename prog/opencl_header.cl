/** @file
 * Inclusion and definition of types and definitions required in the Device code.
 */

//opencl_header.cl

#ifdef cl_amd_printf
// we def out anything using printf if unsupported,
// therefore we also don't need to enable printf
#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

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

