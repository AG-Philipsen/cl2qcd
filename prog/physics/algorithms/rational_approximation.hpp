/** @file
 * Declaration of the Rational_Approximation class
 *
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_RATIONAL_APPROX_
#define _PHYSICS_ALGORITHMS_RATIONAL_APPROX_

#include "../../types.h"
#include"alg_remez.h"
#include <iostream>

namespace physics {
namespace algorithms {

/**
 * Class to calculate all coefficient of a root rational approximation.
 * The function to be approximated is either x^(y/z) or x^(-y/z) and
 * it is understood that the approximation has the form
 * \f[
 *     r(x) = a_0 + sum_{k=1}^{d}  \frac{a_k}{x + b_k}
 * \f]
 * This means that the approximated function must be a ratio of two
 * polynomials of the SAME degree. This is understood in this class.
 * In this way, we can use the Remez algorithm to generate the two
 * polynomials and then calculate the partial fractions expansion.
 * 
 * @note The index of the sum above, that runs from 1 to d, in the code will be
 *       between 0 and d-1 for c++ reasons.
 */
class Rational_Approximation {

public:
	/**
	 * Constructor of the class. It generates the rational approximation calculating
	 * the coefficients a_k and b_k on the basis of the following
	 * @param d The degree of the numerator AND denominator polynomials
	 * @param y The numerator of the exponent
	 * @param z The denominator of the exponent
	 * @param precision The precision that gmp uses
	 * @param low The lower bound of the approximation
	 * @param high The upper bound of the approximation
	 * @param inv If(inv) f_exact=x^(-y/z) else f_exact=x^(y/z)
	 * 
	 * @note The "inv" variable has been set to true by default because in the RHMC
	 *       algorithm we usually have to approximate the function x^(-y/z).
	 * @par
	 * @note The "y" and "z" integer parameter MUST be positive and y%z!=0 to make
	 *       the Remez algorithm work.
	 */
	Rational_Approximation(int d, int y, int z, hmc_float low, hmc_float high,
			        bool inv=true, int precision=100);
	
	/**
	 * Free memory allocated by the constructor
	 */
	~Rational_Approximation();
	
	/*
	 * Rational_Approximation cannot be copied
	 */
	Rational_Approximation& operator=(const Rational_Approximation&) = delete;
	Rational_Approximation(const Rational_Approximation&) = delete;
	Rational_Approximation() = delete;
	
	/**
	 * This function returns the maximum relative error of the approximation: (f_approx-f_exact)/f_exact
	 */
	hmc_float Get_error();
	
	/**
	 * This method returns the order of the approximation, i.e. the degree of the numerator
	 * and denominator polynomials that coincides with the number of partial fractions of the expansion
	 */
	int Get_order();
	
	/**
	 * Method to get the coefficient a0 of the approximation
	 */
	hmc_float Get_a0();
	/**
	 * Method to get the coefficients a of the approximation
	 */
	hmc_float* Get_a();
	/**
	 * Method to get the coefficients b of the approximation
	 */
	hmc_float* Get_b();
	
	/**
	 * Method to get the lower bound of the approximation
	 */
	hmc_float Get_lower_bound();
	/**
	 * Method to get the upper bound of the approximation
	 */
	hmc_float Get_upper_bound();
	
	/**
	 * This method allows the user to print to the shell the information of the approximation
	 */
	friend std::ostream& operator<<(std::ostream&, const Rational_Approximation&); 
	
private:
	bool inv;          /// if(inv) f_exact=x^(-y/z) else f_exact=x^(y/z)
	int d;             /// The degree of the numerator AND denominator polynomials
	int y;             /// The numerator of the exponent
	int z;             /// The denominator of the exponent
	int precision;     /// The precision that gmp uses
	hmc_float low;     /// The lower bound of the approximation
	hmc_float high;    /// The upper bound of the approximation
	hmc_float error;   /// The max relative error of the approximation: (f_approx-f_exact)/f_exact
	/**
	 * The partial fraction expansion takes the form 
	 * r(x) = a0 + sum_{k=0}^{k<d}  a[k] / (x + b[k])
	 */
	hmc_float a0;
	hmc_float *a = new hmc_float[d];
	hmc_float *b = new hmc_float[d];

};


}
}



#endif // _PHYSICS_ALGORITHMS_RATIONAL_APPROX_