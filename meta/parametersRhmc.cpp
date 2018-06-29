/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2017,2018 Alessandro Sciarra
 * Copyright (c) 2018 Francesca Cuteri
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
 */

#include "parametersRhmc.hpp"

using namespace meta;

unsigned int ParametersRhmc::get_md_approx_ord() const noexcept
{
    return molecularDynamicsRationalApproximationOrder;
}
unsigned int ParametersRhmc::get_metro_approx_ord() const noexcept
{
    return metropolisTestRationalApproximationOrder;
}
unsigned int ParametersRhmc::get_findminmax_iteration_block_size() const noexcept
{
    return findMinMaxEigenvalueIterationBlockSize;
}
unsigned int ParametersRhmc::get_findminmax_max() const noexcept
{
    return findMinMaxEigenvalueMaxNumberOfIterations;
}
double ParametersRhmc::get_findminmax_prec() const noexcept
{
    return findMinMaxEigenvaluePrecision;
}
bool ParametersRhmc::get_conservative() const noexcept
{
    return beConservativeInFindMinMaxEigenvalue;
}
double ParametersRhmc::get_num_tastes() const noexcept
{
    return numberOfTastes;
}
unsigned int ParametersRhmc::get_num_tastes_decimal_digits() const noexcept
{
    return numberOfDecimalDigitsInNumberOfTastes;
}
unsigned int ParametersRhmc::get_num_pseudofermions() const noexcept
{
    return numberOfPseudofermions;
}
double ParametersRhmc::get_approx_lower() const noexcept
{
    return lowerBoundForRationalApproximationRange;
}
double ParametersRhmc::get_approx_upper() const noexcept
{
    return upperBoundForRationalApproximationRange;
}
unsigned int ParametersRhmc::get_rhmcsteps() const noexcept
{
    return numberOfRhmcSteps;
}
std::string ParametersRhmc::get_approx_heatbath_file() const noexcept
{
    return heatbathRationalApproximationFilename;
}
std::string ParametersRhmc::get_approx_md_file() const noexcept
{
    return molecularDynamicsRationalApproximationFilename;
}
std::string ParametersRhmc::get_approx_metropolis_file() const noexcept
{
    return metropolisTestRationalApproximationFilename;
}
bool ParametersRhmc::get_read_rational_approximations_from_file() const noexcept
{
    return readRationalApproximationsFromFile;
}

meta::ParametersRhmc::ParametersRhmc() : options("RHMC options")
{
    // clang-format off
	options.add_options()
	("rationalApproxOrderHBAndMD", po::value<unsigned int>(&molecularDynamicsRationalApproximationOrder)->default_value(8), "The order of the rational approximation for the Heat Bath and Molecular Dynamics parts.")
	("rationalApproxOrderMetropolis", po::value<unsigned int>(&metropolisTestRationalApproximationOrder)->default_value(15), "The  order of the rational approximation for the Metropolis part.")
	("findminmaxMaxIteration", po::value<unsigned int>(&findMinMaxEigenvalueMaxNumberOfIterations)->default_value(5000), "The maximum number of iterations in the 'findMinMax' algorithm to find the minimum and the maximum eigenvalues of the fermion matrix operator.")
	("findminmaxResiduumCheckEvery", po::value<unsigned int>(&findMinMaxEigenvalueIterationBlockSize)->default_value(25), "Every how many iteration 'findMinMax' will check the residuum.")
	("findminmaxPrecision", po::value<double>(&findMinMaxEigenvaluePrecision)->default_value(1.e-3), "The precision used in 'findMinMax'.")
	("conservative", po::value<bool>(&beConservativeInFindMinMaxEigenvalue)->default_value(false), "Whether to be conservative in 'findMinMax' (check validity of rational approximation in a wider-than-needed interval). It may affect the correctness of the RHMC, hence use with care.")
	("nTastes", po::value<double>(&numberOfTastes)->default_value(2), "The number of tastes of staggered fermions.")
    ("nTastesDecimalDigits", po::value<unsigned int>(&numberOfDecimalDigitsInNumberOfTastes)->default_value(0), "The number of decimal digits in the number of staggered tastes.")
    ("nPseudoFermions", po::value<unsigned int>(&numberOfPseudofermions)->default_value(1), "The number of pseudo-fermion species in the multiple pseudofermion technique for staggered fermions only.")
	("rationalApproxLowerBound", po::value<double>(&lowerBoundForRationalApproximationRange)->default_value(1.e-5), "The lower bound in the validity interval of the rational approximation.")
	("rationalApproxUpperBound", po::value<double>(&upperBoundForRationalApproximationRange)->default_value(1.), "The upper bound in the validity interval of the rational approximation.")
	("nRhmcSteps", po::value<unsigned int>(&numberOfRhmcSteps)->default_value(10), "The number of RHMC steps (i.e. the number of configuration updates in the Markov chain).")
	("rationalApproxFileHB", po::value<std::string>(&heatbathRationalApproximationFilename)->default_value("Approx_Heatbath"), "The path of the file containing the rational approximation for the Heat Bath part.")
	("rationalApproxFileMD", po::value<std::string>(&molecularDynamicsRationalApproximationFilename)->default_value("Approx_MD"), "The path of the file containing the rational approximation for the Molecular Dynamics part.")
	("rationalApproxFileMetropolis", po::value<std::string>(&metropolisTestRationalApproximationFilename)->default_value("Approx_Metropolis"), "The path of the file containing the rational approximation for the Metropolis part.")
	("readRationalApproxFromFile", po::value<bool>(&readRationalApproximationsFromFile)->default_value(true), "Whether to read rational approximations from file.");
    // clang-format on
}
