/** @file
 * Inclusion and definition of types and definitions required in the Device code.
 */

//opencl_header.cl

#ifdef ENABLE_PRINTF
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

#ifdef _USE_BLOCKED_LOOPS_
#define PARALLEL_FOR(VAR, LIMIT) \
size_t _block_size = (LIMIT + get_global_size(0) - 1) / get_global_size(0); \
for(size_t VAR = get_global_id(0) * _block_size; VAR < (get_global_id(0) + 1) * _block_size && VAR < LIMIT; ++VAR)
#else /* _USE_BLOCKED_LOOPS_ */
#if (NSPACE / 24 ) * 24 == NSPACE /* On Cypress with APP 2.6 we need different strategies for NSPACE % 24 == 0 and NSPACE % 16 == 0 */
#define PARALLEL_FOR(VAR, LIMIT) \
for(size_t VAR = get_global_id(0); VAR < LIMIT; VAR += get_global_size(0))
#else /* NSPACE % 24 == 0 */
#define PARALLEL_FOR(VAR, LIMIT) \
size_t _global_size = get_global_size(0); /* required to avoid performance regression an Cypress with APP 2.6 */ \
for(size_t VAR = get_global_id(0); VAR < LIMIT; VAR += _global_size)
#endif /* NSPACE % 24 == 0 */
#endif /* _USE_BLOCKED_LOOPS_ */
