/** @file
 * Declaration of the Rational_Coefficients and Rational_Approximation classes
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _PHYSICS_ALGORITHMS_RATIONAL_APPROX_
#define _PHYSICS_ALGORITHMS_RATIONAL_APPROX_

#include "../../common_header_files/types.h"
#include"alg_remez.h"
#include "../fermionmatrix/fermionmatrix_stagg.hpp"
#include <iostream>
#include <cmath>
#include "find_minmax_eigenvalue.hpp"
#include "../../host_functionality/logger.hpp"

namespace physics {
namespace algorithms {

/**
 * Class to store all coefficient of a root rational approximation.
 * The approximation has the form
 * \f[
 *     r(x) = a_0 + sum_{k=1}^{d}  \frac{a_k}{x + b_k}
 * \f]
 * so we store here a_0, a_k and b_k. This will be used in the
 * Rational_Approximation class.
 */
  
class Rational_Coefficients {
  
public:
	/**
	 * "Default" constructor of the class. It only allocates memory
	 */
	Rational_Coefficients(const int d);
	
	/**
	 * Constructor of the class. It store the coefficients.
	 *  @param d The degree of approximation, i.e. the size of the vectors.
	 *  @param a0 The first coefficient of the expansion
	 *  @param a The numerator coefficients of the expansion
	 *  @param b The denominator coefficients of the expansion
	 */
	Rational_Coefficients(const int d, const hmc_float a0, const std::vector<hmc_float> a, const std::vector<hmc_float> b);
	
	/**
	 * Default Rational_Coefficients move constructor and move assignment operator
	 * Delete Rational_Coefficients default constructor and copy constructor
	 */
	Rational_Coefficients(Rational_Coefficients&&) = default;
	Rational_Coefficients& operator=(Rational_Coefficients&&) = default;
	Rational_Coefficients& operator=(const Rational_Coefficients&) = delete;
	Rational_Coefficients(const Rational_Coefficients&) = delete;
	Rational_Coefficients() = delete;
    virtual ~Rational_Coefficients(){}

	/**
	 * This method returns the order of the approximation
	 */
    unsigned int Get_order() const;
	/**
	 * Method to get the coefficient a0 of the approximation
	 */
	hmc_float Get_a0() const;
	/**
	 * Method to get the coefficients a of the approximation
	 */
	std::vector<hmc_float> Get_a() const;
	/**
	 * Method to get the coefficients b of the approximation
	 */
	std::vector<hmc_float> Get_b() const;

protected: 
	/**
	 * Method to set the coefficients a0, a and b of the approximation
	 */
	void Set_coeff(const hmc_float v0, const std::vector<hmc_float> v_a, const std::vector<hmc_float> v_b);
  
private:
	int _d;
	/**
	 * The partial fraction expansion takes the form 
	 * r(x) = a0 + sum_{k=0}^{k<d}  a[k] / (x + b[k])
	 */
	hmc_float _a0;
	std::vector<hmc_float> _a;
	std::vector<hmc_float> _b;
  
};


  
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
class Rational_Approximation : public Rational_Coefficients {

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
	 * @par
	 * @note Since this class inherits from Rational_Coefficients and since there are
	 *       no further members that need some memory to be allocated, we do not need
	 *       to implement explictly the destructor (that of Rational_Coefficients is enough)
	 */
	Rational_Approximation(int d, int y, int z, hmc_float low, hmc_float high,
			        bool inv=true, int precision=113);
	
	/**
	 * Constructor of the class from file. It reads the rational approximation parameters
	 * and the coefficients a_0, a_k and b_k.
	 */
	Rational_Approximation(std::string filename);
	
	/*
	 * Rational_Approximation cannot be copied
	 */
	Rational_Approximation& operator=(const Rational_Approximation&) = delete;
	Rational_Approximation(const Rational_Approximation&) = delete;
	Rational_Approximation() = delete;
	
	/**
	 * This function returns the maximum relative error of the approximation: (f_approx-f_exact)/f_exact
	 */
	hmc_float Get_error() const;
	
	/**
	 * Method to get the lower bound of the approximation
	 */
	hmc_float Get_lower_bound() const;
	/**
	 * Method to get the upper bound of the approximation
	 */
	hmc_float Get_upper_bound() const;
	/**
	 * Method to get the exponent of the approximation
	 */
	hmc_float Get_exponent() const;
	
	/**
	 * This method allows to save to text file the rational approximation in a format
	 * compatible to that used to read it in again
	 */
	void Save_rational_approximation(std::string filename);
	
	/**
	 * This method allows the user to print to the shell the information of the approximation
	 */
	friend std::ostream& operator<<(std::ostream&, const Rational_Approximation&); 
	
	/**
	 * This method adapt the coefficients of "this object" to the interval
	 * [lambda_min, lambda_max], where lambda_max and lambda_min are
	 * the maximum and minimum eigenvalues of A (that depends on the gaugefield gf).
	 * This means that the interval where "this object" has been created must be [xmin,1] 
	 * and if this is not the case an error is thrown.
	 * Furthermore, it is checked that xmin is smaller of or equal to lambda_min/lambda_max:
	 * this guarantees the output of this method to be reliable. If this check fails
	 * another error is thrown.
	 */
	Rational_Coefficients Rescale_Coefficients(const hmc_float minEigenvalue, const hmc_float maxEigenvalue) const;

private:
	bool inv;          /// if(inv) f_exact=x^(-y/z) else f_exact=x^(y/z)
	int y;             /// The numerator of the exponent
	int z;             /// The denominator of the exponent
	int precision;     /// The precision that gmp uses
	hmc_float low;     /// The lower bound of the approximation
	hmc_float high;    /// The upper bound of the approximation
	hmc_float error;   /// The max relative error of the approximation: (f_approx-f_exact)/f_exact

};


}
}

#endif // _PHYSICS_ALGORITHMS_RATIONAL_APPROX_
