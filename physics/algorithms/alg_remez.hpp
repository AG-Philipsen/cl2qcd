/**
 * Copyright (c) 2013,2018 Alessandro Sciarra
 * Copyright (c) 2013 Matthias Bach
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 *
 * NOTE: The code contained in this file was developed by external developers
 *       and the copyright and license statements above refer to the work
 *       that was done to include the third party code into CL2QCD.
 */

/** @file
 *
 * Copyright (c) Mike Clark - 25th May 2005
 *
 * alg_remez.h
 *
 * AlgRemez is an implementation of the Remez algorithm, which in this
 * case is used for generating the optimal nth root rational
 * approximation.
 *
 * Note this class requires the gnu multiprecision (GNU MP) library.
 *
 * @TODO Add documentation with doxygen (so far see code)
 */

#ifndef _PHYSICS_ALGORITHMS_ALG_REMEZ_
#define _PHYSICS_ALGORITHMS_ALG_REMEZ_

#include "bigfloat.hpp"

#define JMAX 10000 //Maximum number of iterations of Newton's approximation
#define SUM_MAX 10 // Maximum number of terms in exponential

namespace physics {
namespace algorithms {

class AlgRemez
{
 private:
  char *cname;

  // The approximation parameters
  bigfloat *param, *roots, *poles;
  bigfloat norm;

  // The numerator and denominator degree (n=d)
  int n, d;

  // The bounds of the approximation
  bigfloat apstrt, apwidt, apend;

  // the numerator and denominator of the power we are approximating
  unsigned long power_num;
  unsigned long power_den;

  // Flag to determine whether the arrays have been allocated
  int alloc;

  // Flag to determine whether the roots have been found
  int foundRoots;

  // Variables used to calculate the approximation
  int nd1, iter;
  bigfloat *xx, *mm, *step;
  bigfloat delta, spread, tolerance;

  // The exponential summation coefficients
  bigfloat *a;
  int *a_power;
  int a_length;

  // The number of equations we must solve at each iteration (n+d+1)
  int neq;

  // The precision of the GNU MP library
  long prec;

  // Initial values of maximal and minmal errors
  void initialGuess();

  // Solve the equations
  void equations();

  // Search for error maxima and minima
  void search(bigfloat *step);

  // Initialise step sizes
  void stpini(bigfloat *step);

  // Calculate the roots of the approximation
  int root();

  // Evaluate the polynomial
  bigfloat polyEval(bigfloat x, bigfloat *poly, long size);
  //complex_bf polyEval(complex_bf x, complex_bf *poly, long size);

  // Evaluate the differential of the polynomial
  bigfloat polyDiff(bigfloat x, bigfloat *poly, long size);
  //complex_bf polyDiff(complex_bf x, complex_bf *poly, long size);

  // Newton's method to calculate roots
  bigfloat rtnewt(bigfloat *poly, long i, bigfloat x1, bigfloat x2, bigfloat xacc);
  //complex_bf rtnewt(complex_bf *poly, long i, bigfloat x1, bigfloat x2, bigfloat xacc);

  // Evaluate the partial fraction expansion of the rational function
  // with res roots and poles poles.  Result is overwritten on input
  // arrays.
  void pfe(bigfloat *res, bigfloat* poles, bigfloat norm);

  // Calculate function required for the approximation
  bigfloat func(bigfloat x);

  // Compute size and sign of the approximation error at x
  bigfloat getErr(bigfloat x, int *sign);

  // Solve the system AX=B
  int simq(bigfloat *A, bigfloat *B, bigfloat *X, int n);

  // Free memory and reallocate as necessary
  void allocate(int num_degree, int den_degree);

  // Evaluate the rational form P(x)/Q(x) using coefficients from the
  // solution vector param
  bigfloat approx(bigfloat x);

 public:

  // Constructor
  AlgRemez(double lower, double upper, long prec);

  // Destructor
  virtual ~AlgRemez();

  // Reset the bounds of the approximation
  void setBounds(double lower, double upper);

  // Generate the rational approximation x^(pnum/pden)
  double generateApprox(int num_degree, int den_degree,
			unsigned long power_num, unsigned long power_den,
			int a_len, double* a_param, int* a_pow);
  double generateApprox(int num_degree, int den_degree,
			unsigned long power_num, unsigned long power_den);
  double generateApprox(int degree, unsigned long power_num,
			unsigned long power_den);

  // Return the partial fraction expansion of the approximation x^(pnum/pden)
  int getPFE(double *res, double *pole, double *norm);

  // Return the partial fraction expansion of the approximation x^(-pnum/pden)
  int getIPFE(double *res, double *pole, double *norm);

  // Evaluate the rational form P(x)/Q(x) using coefficients from the
  // solution vector param
  double evaluateApprox(double x);

  // Evaluate the rational form Q(x)/P(x) using coefficients from the
  // solution vector param
  double evaluateInverseApprox(double x);

  // Calculate function required for the approximation
  double evaluateFunc(double x);

  // Calculate inverse function required for the approximation
  double evaluateInverseFunc(double x);

};

}
}
#endif  // _PHYSICS_ALGORITHMS_ALG_REMEZ_
