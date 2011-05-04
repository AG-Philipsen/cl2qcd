/** @file
 * Operations on 3x3 matrices
 */
#ifndef _MATRIXH_
#define _MATRIXH_

#include "types.h"
#include "hmcerrs.h"
#include "host_operations_complex.h"

/**
 * Multiplies two 3x3 matrices
 * @param[out] out The matrix into which to store the multiplication result
 * @param[in] p Left matrix for the multiplication
 * @param[in] q Right matrix for the multiplication
 * @return Error code as defined in hmcerrs.h
 */
hmc_error multiply_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);
/**
 * Adds two 3x3 matrices
 * @param[out] out The matrix into which to store the sum result
 * @param[in] p Left matrix for the summation
 * @param[in] q Right matrix for the summation
 * @return Error code as defined in hmcerrs.h
 */
hmc_error add_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);
/**
 * Subtracts two 3x3 matrices
 * @param[out] out The matrix into which to store the difference result
 * @param[in] p Left matrix for the subtraction
 * @param[in] q Right matrix for the subtraction
 * @return Error code as defined in hmcerrs.h
 */
hmc_error substract_3x3matrix (hmc_3x3matrix *out, hmc_3x3matrix *p, hmc_3x3matrix *q);

#endif