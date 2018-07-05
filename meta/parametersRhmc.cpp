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
unsigned int ParametersRhmc::get_rhmcsteps() const noexcept
{
    return numberOfRhmcSteps;
}

meta::ParametersRhmc::ParametersRhmc()
    : numberOfTastes(2)
    , numberOfDecimalDigitsInNumberOfTastes(0)
    , numberOfPseudofermions(1)
    , numberOfRhmcSteps(10)
    , options("RHMC options")
{
    // clang-format off
    options.add_options()
    ("nTastes", po::value<double>(&numberOfTastes)->default_value(numberOfTastes), "The number of tastes of staggered fermions.")
    ("nTastesDecimalDigits", po::value<unsigned int>(&numberOfDecimalDigitsInNumberOfTastes)->default_value(numberOfDecimalDigitsInNumberOfTastes), "The number of decimal digits in the number of staggered tastes.")
    ("nPseudoFermions", po::value<unsigned int>(&numberOfPseudofermions)->default_value(numberOfPseudofermions), "The number of pseudo-fermion species in the multiple pseudofermion technique for staggered fermions only.")
    ("nRhmcSteps", po::value<unsigned int>(&numberOfRhmcSteps)->default_value(numberOfRhmcSteps), "The number of RHMC steps (i.e. the number of configuration updates in the Markov chain).");
    // clang-format on
}
