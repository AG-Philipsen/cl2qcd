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

//Rational_Coefficients class

physics::algorithms::Rational_Coefficients::Rational_Coefficients(const int d) : _d(d), _a(new hmc_float[d]),
_b(new hmc_float[d]) { }

physics::algorithms::Rational_Coefficients::Rational_Coefficients(const int d, const hmc_float a0, const hmc_float* a, const hmc_float* b) : _d(d), _a0(a0)
{
	//Allocate and initialize
	_a = new hmc_float[d];
	_b = new hmc_float[d];
	for(int i=0; i<d; i++){
		_a[i] = a[i];
		_b[i] = b[i];
	}
}

physics::algorithms::Rational_Coefficients::~Rational_Coefficients()
{
	delete[] _a;
	delete[] _b;
}

int physics::algorithms::Rational_Coefficients::Get_order() const
{
	return _d;
}


hmc_float physics::algorithms::Rational_Coefficients::Get_a0() const
{
	return _a0;
}


hmc_float* physics::algorithms::Rational_Coefficients::Get_a() const
{
	return _a;
}


hmc_float* physics::algorithms::Rational_Coefficients::Get_b() const
{
	return _b;
}

void physics::algorithms::Rational_Coefficients::Set_a0(hmc_float v) 
{
	_a0 = v;
}


void physics::algorithms::Rational_Coefficients::Set_a(hmc_float* v) 
{
	for(int i=0; i<_d; i++)
		_a[i] = v[i];
}


void physics::algorithms::Rational_Coefficients::Set_b(hmc_float* v) 
{
	for(int i=0; i<_d; i++)
		_b[i] = v[i];
}


//Rational_Approximation class

physics::algorithms::Rational_Approximation::Rational_Approximation(int d, int y, int z, hmc_float low, hmc_float high, bool inv, int precision) : Rational_Coefficients(d), inv(inv), y(y), z(z), precision(precision), low(low), high(high)
{
	//Checks on the exponent numerator and denominator
	if(y<=0 || z<=0 || y%z==0)
	  throw std::invalid_argument("Exponent of the rational approximation not allowed!");
	if(inv==false){
	  logger.info() << "Calculating the rational approx. of degree " << 
	                    d << " of the function x^("<< y << "/" << z << ")...";
	}else{
	  logger.info() << "Calculating the rational approx. of degree " << 
	                    d << " of the function x^(-"<< y << "/" << z << ")...";
	}
	// Instantiate the Remez class
	physics::algorithms::AlgRemez remez(low, high, precision);
	// Generate the required approximation into the form P(x)/Q(x)
	// Remark that the method generateApprox returns the maximum relative error occurred
	// and it is the same to approximate both x^{y/z}  and x^{-y/z}
	error = remez.generateApprox(d,d,y,z);
	// Find the partial fraction expansion of the approximation either
	// to the function x^{y/z} or to the function x^{-y/z}
	hmc_float *a_tmp = new hmc_float[d];
	hmc_float *b_tmp = new hmc_float[d];
	hmc_float *a0_tmp = new hmc_float;
	if(inv==false)
	  remez.getPFE(a_tmp, b_tmp, a0_tmp);
	else
	  remez.getIPFE(a_tmp, b_tmp, a0_tmp);
	
	Set_a(a_tmp);
	Set_b(b_tmp);
	Set_a0(*a0_tmp);
	
	delete[] a_tmp;
	delete[] b_tmp;
	delete a0_tmp;
}


hmc_float physics::algorithms::Rational_Approximation::Get_error()
{
	return error;
}


hmc_float physics::algorithms::Rational_Approximation::Get_lower_bound()
{
	return low;
}

hmc_float physics::algorithms::Rational_Approximation::Get_upper_bound()
{
	return high;
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
	os << "\t\t ++\t         Approximation order = " << approx.Get_order() << "\t\t\t\t++\n";
	os << "\t\t ++\t  Range of the approximation = [" << approx.low << "," << approx.high <<"]\t\t\t++\n";
	os << "\t\t ++\t Precision of the arithmetic = " << approx.precision << "\t\t\t\t++\n";
	os << "\t\t ++\t      Maximum relative error = " << approx.error << "\t\t++\n";
	os.precision(16);
	os.setf( std::ios::fixed, std:: ios::floatfield);
	os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
	os << "\t\t ++        a0 = "<< approx.Get_a0();
	os << "\t\t\t\t\t++\n";
	for(int i=0; i<approx.Get_order(); i++){
	  if(i<9){
	    os << "\t\t ++      a[" << i+1 << "] = "<< approx.Get_a()[i] << "       ";
	    os << " b[" << i+1 << "] = "<< approx.Get_b()[i]<<"\t++\n";
	  }else{
	    os << "\t\t ++     a[" << i+1 << "] = "<< approx.Get_a()[i] << "       ";
	    os << "b[" << i+1 << "] = "<< approx.Get_b()[i]<<"\t++\n";
	  }
	}
	os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
	os << "\t\t +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
	return os;
}
}
}



