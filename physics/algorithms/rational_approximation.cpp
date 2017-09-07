/** @file
 * Implementation of the Rational_Approximation class
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

#include "rational_approximation.hpp"
//#include "alg_remez.h"
//#include "../../types.h"
//#include "../../exceptions.h"
//#include "../../logger.hpp"
//#include <fstream>
//#include <cmath>

//Rational_Coefficients class

physics::algorithms::Rational_Coefficients::Rational_Coefficients(const int d)
        : _d(d), _a0()
{
    //This constructor can appear unnecessary. In fact, reserving memory we avoid
    //reallocations pushing back elements in the vectors, since in this way we are
    //sure that the capacity of the vectors _a and _b is at least d.
    _a.reserve(d);
    _b.reserve(d);
}

physics::algorithms::Rational_Coefficients::Rational_Coefficients(const int d, const hmc_float a0, const std::vector<hmc_float> a,
        const std::vector<hmc_float> b)
        : _d(d), _a0(a0)
{
    //Reserve and initialize
    _a.reserve(d);
    _b.reserve(d);
    for (int i = 0; i < d; i++) {
        _a.push_back(a[i]);
        _b.push_back(b[i]);
    }
}

unsigned int physics::algorithms::Rational_Coefficients::Get_order() const
{
	return _d;
}

hmc_float physics::algorithms::Rational_Coefficients::Get_a0() const
{
    return _a0;
}

std::vector<hmc_float> physics::algorithms::Rational_Coefficients::Get_a() const
{
    return _a;
}

std::vector<hmc_float> physics::algorithms::Rational_Coefficients::Get_b() const
{
    return _b;
}

void physics::algorithms::Rational_Coefficients::Set_coeff(const hmc_float v0, const std::vector<hmc_float> v_a, const std::vector<hmc_float> v_b)
{
    if(v_a.size() != v_b.size())
        throw std::invalid_argument("Vectors with different sizes passed to Rational_Coefficients::Set_coeff!");
    if(_d != (int) v_a.size()) {
        logger.debug() << "Rational_Coefficients::Set_coeff changed the order of the instance.";
        _d = v_a.size();
    }
    _a0 = v0;
    _a = v_a;   //Recall that the assignment operator between vectors "Assigns new contents to the
    _b = v_b;   //container, replacing its current contents, and modifying its size accordingly".
}

//Rational_Approximation class

physics::algorithms::Rational_Approximation::Rational_Approximation(int d, int y, int z, hmc_float low, hmc_float high, bool inv, int precision)
        : Rational_Coefficients(d), inv(inv), y(y), z(z), precision(precision), low(low), high(high)
{
    //Checks on the exponent numerator and denominator
    if(y <= 0 || z <= 0 || y % z == 0)
        throw std::invalid_argument("Exponent of the rational approximation not allowed!");
    if(inv == false) {
        logger.info() << "Calculating the rational approx. of degree " << d << " of the function x^(" << y << "/" << z << ")...";
    } else {
        logger.info() << "Calculating the rational approx. of degree " << d << " of the function x^(-" << y << "/" << z << ")...";
    }
    // Instantiate the Remez class
    physics::algorithms::AlgRemez remez(low, high, precision);
    // Generate the required approximation into the form P(x)/Q(x)
    // Remark that the method generateApprox returns the maximum relative error occurred
    // and it is the same to approximate both x^{y/z}  and x^{-y/z}
    error = remez.generateApprox(d, d, y, z);
    // Find the partial fraction expansion of the approximation either
    // to the function x^{y/z} or to the function x^{-y/z}
    double *a_tmp = new hmc_float[d];
    double *b_tmp = new hmc_float[d];
    double *a0_tmp = new hmc_float;
    if(inv == false)
        remez.getPFE(a_tmp, b_tmp, a0_tmp);
    else
        remez.getIPFE(a_tmp, b_tmp, a0_tmp);

    //I need std::vector to use Set_coeff
    std::vector<hmc_float> av_tmp(a_tmp, a_tmp + d);
    std::vector<hmc_float> bv_tmp(b_tmp, b_tmp + d);

    Set_coeff(*a0_tmp, av_tmp, bv_tmp);

    delete[] a_tmp;
    delete[] b_tmp;
    delete a0_tmp;
}

physics::algorithms::Rational_Approximation::Rational_Approximation(std::string filename)
        : Rational_Coefficients(0), inv(), y(), z(), precision(), low(), high(), error()
{
    std::fstream file;
    file.open(filename.c_str());
    if(!file) {
        logger.fatal() << "Unable to open the file with Rational_Approximation data!";
        throw File_Exception(filename);
    }
    int d;
    hmc_float a0, tmp;
    std::vector<hmc_float> a, b;
    file >> d >> y >> z >> low >> high >> inv >> precision >> error >> a0;
    for (int i = 0; i < d; i++) {
        file >> tmp;
        a.push_back(tmp);
        file >> tmp;
        b.push_back(tmp);
    }
    file >> tmp;
    if(!(file.peek() == EOF && file.eof()))
        logger.warn() << "The file with Rational_Approximation data contains additional stuff! Check it!";
    file.close();
    Set_coeff(a0, a, b);
}

void physics::algorithms::Rational_Approximation::Save_rational_approximation(std::string filename)
{
    std::fstream outputToFile;
    outputToFile.open(filename.c_str(), std::ios::out);
    if(!outputToFile.is_open())
        throw File_Exception(filename);
    outputToFile.precision(15);
    outputToFile << Get_order() << "     " << y << " " << z << "     " << low << " " << high << "     ";
    outputToFile << inv << " " << precision << "     ";
    outputToFile.setf(std::ios::scientific, std::ios::floatfield);
    outputToFile << error << "\n\n" << Get_a0() << "\n";
    std::vector<hmc_float> a_tmp = Get_a();
    std::vector<hmc_float> b_tmp = Get_b();
    for(unsigned int i = 0; i < Get_order(); i++)
        outputToFile << a_tmp[i] << "\t\t" << b_tmp[i] << "\n";
    outputToFile.close();
}

hmc_float physics::algorithms::Rational_Approximation::Get_error() const
{
    return error;
}

hmc_float physics::algorithms::Rational_Approximation::Get_lower_bound() const
{
    return low;
}

hmc_float physics::algorithms::Rational_Approximation::Get_upper_bound() const
{
    return high;
}

hmc_float physics::algorithms::Rational_Approximation::Get_exponent() const
{
    if(inv)
        return -1 * ((hmc_float) y) / z;
    else
        return ((hmc_float) y) / z;
}

physics::algorithms::Rational_Coefficients physics::algorithms::Rational_Approximation::Rescale_Coefficients(const hmc_float minEigenvalue, const hmc_float maxEigenvalue) const
{
	if(high != 1)
		throw std::invalid_argument("Upper bound different from 1 in rescale_coefficients!");

	if(low > minEigenvalue/maxEigenvalue)
	        throw Print_Error_Message("Rational_Approximation does not respect lower_bound <= min/max");

	int ord = Get_order();
	hmc_float exp = Get_exponent();
	hmc_float a0_new = Get_a0();
	std::vector<hmc_float> a_new = Get_a();
	std::vector<hmc_float> b_new = Get_b();

	hmc_float tmp = pow(maxEigenvalue, exp);
	a0_new *= tmp;
	for(int i=0; i<ord; i++){
		a_new[i] *= (tmp*maxEigenvalue);
		b_new[i] *= maxEigenvalue;
	}

	Rational_Coefficients out(ord, a0_new, a_new, b_new);
	return out;
}

namespace physics {
    namespace algorithms {

        std::ostream& operator <<(std::ostream &os, const Rational_Approximation &approx)
        {
            if(approx.inv) {
                os << "\n\t\t +++++++++++++++++++ Approximation to f(x) = x^(-" << approx.y << "/" << approx.z;
                os << ") ++++++++++++++++++++\n";
            } else {
                os << "\n\t\t ++++++++++++++++++++ Approximation to f(x) = x^(" << approx.y << "/" << approx.z;
                os << ") ++++++++++++++++++++\n";
            }

            os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
            os << "\t\t ++\t         Approximation order = " << approx.Get_order() << "\t\t\t\t++\n";
            os << "\t\t ++\t  Range of the approximation = [" << approx.low << "," << approx.high << "]\t\t\t++\n";
            os << "\t\t ++\t Precision of the arithmetic = " << approx.precision << "\t\t\t\t++\n";
            os << "\t\t ++\t      Maximum relative error = " << approx.error << "\t\t++\n";
            os.precision(16);
            os.setf(std::ios::fixed, std::ios::floatfield);
            os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
            os << "\t\t ++        a0 = " << approx.Get_a0();
            os << "\t\t\t\t\t++\n";
            for(unsigned int i = 0; i < approx.Get_order(); i++) {
                if(i < 9) {
                    os << "\t\t ++      a[" << i + 1 << "] = " << approx.Get_a()[i] << "       ";
                    os << " b[" << i + 1 << "] = " << approx.Get_b()[i] << "\t++\n";
                } else {
                    os << "\t\t ++     a[" << i + 1 << "] = " << approx.Get_a()[i] << "       ";
                    os << "b[" << i + 1 << "] = " << approx.Get_b()[i] << "\t++\n";
                }
            }
            os << "\t\t ++\t\t\t\t\t\t\t\t\t++\n";
            os << "\t\t +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
            return os;
        }
    }
}

