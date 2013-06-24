/** @file
 * Tests of the Rational_Approximation algorithm
 */

#include "rational_approximation.hpp"
#include "../../logger.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::rational_approximation
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::algorithms;
	
	Rational_Approximation approx(6,1,2,1e-5,1);
	logger.info() << approx;
}

#if 0
BOOST_AUTO_TEST_CASE(rational_approximation)
{
  
  int i=0;
  int n=6; // The degree of the numerator polynomial
  int d=6; // The degree of the denominator polynomial
  int y=1; // The numerator of the exponent
  int z=2; // The denominator of the exponent
  int precision=40; // The precision that gmp uses
  double lambda_low=0.0004;
  double lambda_high=64; // The bounds of the approximation

  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res = new double[n];
  double *pole = new double[d];

  double bulk = exp(0.5*(log(lambda_low)+log(lambda_high)));
  char FORCE_FILE[100], ENERGY_FILE[100];
  sprintf(FORCE_FILE, "force_%d_%d_%d_%d_%f.dat", y, z, d, n, bulk);
  sprintf(ENERGY_FILE, "energy_%d_%d_%d_%d_%f.dat", y, z, d, n, bulk);

  // Instantiate the Remez class
  physics::algorithms::AlgRemez remez(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez.generateApprox(n,d,y,z);

  //FILE *output = fopen("approx.dat", "w");

  //fprintf(output, "Approximation to f(x) = x^(%d/%d)\n\n", y, z);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  remez.getPFE(res,pole,&norm);
  
  /*
  fprintf(output, "alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	    i+1, res[i], i+1, pole[i]);
  }
  */
  
  // Find pfe of inverse function
  remez.getIPFE(res,pole,&norm);
//   fprintf(output, "\nApproximation to f(x) = x^(-%d/%d)\n\n", y, z);
//   fprintf(output, "alpha[0] = %18.16e\n", norm);
  printf("\nApproximation to f(x) = x^(-%d/%d)\n\n", y, z);
  printf("alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
//     fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
// 	    i+1, res[i], i+1, pole[i]);
     printf("alpha[%d] = %18.16e, beta[%d] = %18.16e\n", i+1, res[i], i+1, pole[i]);
  }

  //fclose(output);

  /*
  FILE *error_file = fopen("error.dat", "w");
  for (double x=lambda_low; x<lambda_high; x*=1.01) {
    double f = remez.evaluateFunc(x);
    double r = remez.evaluateApprox(x);
    fprintf(error_file,"%e %e\n", x,  (r - f)/f);
  }
  fclose(error_file);
  */

  delete res;
  delete pole;
}     
#endif