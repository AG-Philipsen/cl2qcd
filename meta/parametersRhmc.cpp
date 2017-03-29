/** @file
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

meta::ParametersRhmc::ParametersRhmc()
	: options("RHMC options")
{
	options.add_options()
	("md_approx_ord", po::value<unsigned int>(&molecularDynamicsRationalApproximationOrder)->default_value(8))
	("metro_approx_ord", po::value<unsigned int>(&metropolisTestRationalApproximationOrder)->default_value(15))
	("findminmax_max", po::value<unsigned int>(&findMinMaxEigenvalueMaxNumberOfIterations)->default_value(5000))
	("findminmax_iteration_block_size", po::value<unsigned int>(&findMinMaxEigenvalueIterationBlockSize)->default_value(25), "find_minmax will check the residuum only every N iterations")
	("findminmax_prec", po::value<double>(&findMinMaxEigenvaluePrecision)->default_value(1.e-3))
	("conservative", po::value<bool>(&beConservativeInFindMinMaxEigenvalue)->default_value(false))
	("num_tastes", po::value<double>(&numberOfTastes)->default_value(2))
    ("num_tastes_decimal_digits", po::value<unsigned int>(&numberOfDecimalDigitsInNumberOfTastes)->default_value(0))
    ("num_pseudofermions", po::value<unsigned int>(&numberOfPseudofermions)->default_value(1), "ONLY for staggered fermions!")
	("approx_lower", po::value<double>(&lowerBoundForRationalApproximationRange)->default_value(1.e-5))
	("approx_upper", po::value<double>(&upperBoundForRationalApproximationRange)->default_value(1.))
	("rhmcsteps", po::value<unsigned int>(&numberOfRhmcSteps)->default_value(10))
	("approx_heatbath_file", po::value<std::string>(&heatbathRationalApproximationFilename)->default_value("Approx_Heatbath"))
	("approx_md_file", po::value<std::string>(&molecularDynamicsRationalApproximationFilename)->default_value("Approx_MD"))
	("approx_metropolis_file", po::value<std::string>(&metropolisTestRationalApproximationFilename)->default_value("Approx_Metropolis"))
	("read_rational_approximations_from_file", po::value<bool>(&readRationalApproximationsFromFile)->default_value(true));
}

meta::ParametersRhmc::~ParametersRhmc() = default;

po::options_description & meta::ParametersRhmc::getOptions()
{
	return options;
}
