/** @file
 *
 * Copyright (c) 2018 Alessandro Sciarra
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

#include "parametersRationalApproximation.hpp"

using namespace meta;

unsigned int ParametersRationalApproximation::get_md_approx_ord() const noexcept
{
    return molecularDynamicsRationalApproximationOrder;
}
unsigned int ParametersRationalApproximation::get_metro_approx_ord() const noexcept
{
    return metropolisTestRationalApproximationOrder;
}
unsigned int ParametersRationalApproximation::get_findminmax_iteration_block_size() const noexcept
{
    return findMinMaxEigenvalueIterationBlockSize;
}
unsigned int ParametersRationalApproximation::get_findminmax_max() const noexcept
{
    return findMinMaxEigenvalueMaxNumberOfIterations;
}
double ParametersRationalApproximation::get_findminmax_prec() const noexcept
{
    return findMinMaxEigenvaluePrecision;
}
bool ParametersRationalApproximation::get_conservative() const noexcept
{
    return beConservativeInFindMinMaxEigenvalue;
}
double ParametersRationalApproximation::get_approx_lower() const noexcept
{
    return lowerBoundForRationalApproximationRange;
}
double ParametersRationalApproximation::get_approx_upper() const noexcept
{
    return upperBoundForRationalApproximationRange;
}
std::string ParametersRationalApproximation::get_approx_heatbath_file() const noexcept
{
    return heatbathRationalApproximationFilename;
}
std::string ParametersRationalApproximation::get_approx_md_file() const noexcept
{
    return molecularDynamicsRationalApproximationFilename;
}
std::string ParametersRationalApproximation::get_approx_metropolis_file() const noexcept
{
    return metropolisTestRationalApproximationFilename;
}
bool ParametersRationalApproximation::get_read_rational_approximations_from_file() const noexcept
{
    return readRationalApproximationsFromFile;
}

meta::ParametersRationalApproximation::ParametersRationalApproximation()
    : molecularDynamicsRationalApproximationOrder(8)
    , metropolisTestRationalApproximationOrder(15)
    , findMinMaxEigenvalueIterationBlockSize(25)
    , findMinMaxEigenvalueMaxNumberOfIterations(5000)
    , findMinMaxEigenvaluePrecision(1.e-3)
    , beConservativeInFindMinMaxEigenvalue(false)
    , lowerBoundForRationalApproximationRange(1.e-5)
    , upperBoundForRationalApproximationRange(1.)
    , readRationalApproximationsFromFile(true)
    , heatbathRationalApproximationFilename("Approx_Heatbath")
    , molecularDynamicsRationalApproximationFilename("Approx_MD")
    , metropolisTestRationalApproximationFilename("Approx_Metropolis")
    , options("Rational Approximation options")
{
    // clang-format off
    options.add_options()
    ("rationalApproxOrderHBAndMD", po::value<unsigned int>(&molecularDynamicsRationalApproximationOrder)->default_value(molecularDynamicsRationalApproximationOrder), "The order of the rational approximation for the Heat Bath and Molecular Dynamics parts.")
    ("rationalApproxOrderMetropolis", po::value<unsigned int>(&metropolisTestRationalApproximationOrder)->default_value(metropolisTestRationalApproximationOrder), "The  order of the rational approximation for the Metropolis part.")
    ("findminmaxMaxIteration", po::value<unsigned int>(&findMinMaxEigenvalueMaxNumberOfIterations)->default_value(findMinMaxEigenvalueMaxNumberOfIterations), "The maximum number of iterations in the 'findMinMax' algorithm to find the minimum and the maximum eigenvalues of the fermion matrix operator.")
    ("findminmaxResiduumCheckEvery", po::value<unsigned int>(&findMinMaxEigenvalueIterationBlockSize)->default_value(findMinMaxEigenvalueIterationBlockSize), "Every how many iteration 'findMinMax' will check the residuum.")
    ("findminmaxPrecision", po::value<double>(&findMinMaxEigenvaluePrecision)->default_value(findMinMaxEigenvaluePrecision), "The precision used in 'findMinMax'.")
    ("conservative", po::value<bool>(&beConservativeInFindMinMaxEigenvalue)->default_value(beConservativeInFindMinMaxEigenvalue), "Whether to be conservative in 'findMinMax' (check validity of rational approximation in a wider-than-needed interval). It may affect the correctness of the RHMC, hence use with care.")
    ("rationalApproxLowerBound", po::value<double>(&lowerBoundForRationalApproximationRange)->default_value(lowerBoundForRationalApproximationRange), "The lower bound in the validity interval of the rational approximation.")
    ("rationalApproxUpperBound", po::value<double>(&upperBoundForRationalApproximationRange)->default_value(upperBoundForRationalApproximationRange), "The upper bound in the validity interval of the rational approximation.")
    ("rationalApproxFileHB", po::value<std::string>(&heatbathRationalApproximationFilename)->default_value(heatbathRationalApproximationFilename), "The path of the file containing the rational approximation for the Heat Bath part.")
    ("rationalApproxFileMD", po::value<std::string>(&molecularDynamicsRationalApproximationFilename)->default_value(molecularDynamicsRationalApproximationFilename), "The path of the file containing the rational approximation for the Molecular Dynamics part.")
    ("rationalApproxFileMetropolis", po::value<std::string>(&metropolisTestRationalApproximationFilename)->default_value(metropolisTestRationalApproximationFilename), "The path of the file containing the rational approximation for the Metropolis part.")
    ("readRationalApproxFromFile", po::value<bool>(&readRationalApproximationsFromFile)->default_value(readRationalApproximationsFromFile), "Whether to read rational approximations from file.");
    // clang-format on
}
