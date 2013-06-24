/** @file
 * Implementation of the Rational_Approximation class
 *
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "rational_approximation.hpp"
#include"alg_remez.h"
#include "../../types.h"
#include "../../exceptions.h"
#include "../../logger.hpp"

physics::algorithms::Rational_Approximation::Rational_Approximation(int d, int y, int z, hmc_float low, hmc_float high, bool inv, int precision) : inv(inv), d(d), y(y), z(z), precision(precision), low(low), high(high)
{
	//Checks on the exponent numerator and denominator
	if(y<=0 || z<=0 || y%z==0)
	  throw std::invalid_argument("Exponent of the rational approximation not allowed!");
	if(inv==false)
	  logger.info() << "Calculating the rational approximation of the function x^("<< y << "/" << z << ")...";
	else
	  logger.info() << "Calculating the rational approximation of the function x^(-"<< y << "/" << z << ")...";
	//Allocate a and b
	a = new hmc_float[d];
	b = new hmc_float[d];
	// Instantiate the Remez class
	physics::algorithms::AlgRemez remez(low, high, precision);
	// Generate the required approximation into the form P(x)/Q(x)
	// Remark that the method generateApprox returns the maximum relative error occurred
	// and it is the same to approximate both x^{y/z}  and x^{-y/z}
	error = remez.generateApprox(d,d,y,z);
	// Find the partial fraction expansion of the approximation either
	// to the function x^{y/z} or to the function x^{-y/z}
	if(inv==false)
	  remez.getPFE(a,b,&a0);
	else
	  remez.getIPFE(a,b,&a0);
}


physics::algorithms::Rational_Approximation::~Rational_Approximation()
{
	delete[] a;
	delete[] b;
}


hmc_float physics::algorithms::Rational_Approximation::Get_error()
{
	return error;
}


int physics::algorithms::Rational_Approximation::Get_order()
{
	return d;
}


hmc_float physics::algorithms::Rational_Approximation::Get_a0()
{
	return a0;
}


hmc_float* physics::algorithms::Rational_Approximation::Get_a()
{
	return a;
}


hmc_float* physics::algorithms::Rational_Approximation::Get_b()
{
	return b;
}


namespace physics {
namespace algorithms {
  
std::ostream& operator <<(std::ostream &os, const Rational_Approximation &approx)
{
	if(approx.inv){
	  os << "\n\t\t +++++++++++++++++++ Approximation to f(x) = x^(-" << approx.y << "/" << approx.z;
	  os << ") ++++++++++++++++++++\n";
	}else{
	  os << "\n\t\t ++++++++++++++++++++ Approximation to f(x) = x^(" << approx.y << "/" << approx.z;
	  os << ") ++++++++++++++++++++\n";
	}
	os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
	os << "\t\t ++\t         Approximation order = " << approx.d << "\t\t\t\t++\n";
	os << "\t\t ++\t  Range of the approximation = [" << approx.low << "," << approx.high <<"]\t\t\t++\n";
	os << "\t\t ++\t Precision of the arithmetic = " << approx.precision << "\t\t\t\t++\n";
	os << "\t\t ++\t      Maximum relative error = " << approx.error << "\t\t\t++\n";
	os.precision(16);
	os.setf( std::ios::fixed, std:: ios::floatfield);
	os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
	os << "\t\t ++        a0 = "<< approx.a0;
	os << "\t\t\t\t\t++\n";
	for(int i=0; i<approx.d; i++){
	  if(i<9){
	    os << "\t\t ++      a[" << i+1 << "] = "<< approx.a[i] << "       ";
	    os << " b[" << i+1 << "] = "<< approx.b[i]<<"\t++\n";
	  }else{
	    os << "\t\t ++     a[" << i+1 << "] = "<< approx.a[i] << "       ";
	    os << "b[" << i+1 << "] = "<< approx.b[i]<<"\t++\n";
	  }
	}
	os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
	os << "\t\t +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
	return os;
}
}
}



