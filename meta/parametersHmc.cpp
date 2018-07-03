/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
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

#include "parametersHmc.hpp"

#include "../executables/exceptions.hpp"
#include "../host_functionality/logger.hpp"

#include <boost/algorithm/string.hpp>
#include <stdexcept>

using namespace meta;

int ParametersHmc::get_hmcsteps() const noexcept
{
    return hmcsteps;
}

bool ParametersHmc::get_use_gauge_only() const noexcept
{
    return use_gauge_only;
}

bool ParametersHmc::get_use_mp() const noexcept
{
    return use_mp;
}

ParametersHmc::ParametersHmc() : options("HMC options")
{
    // clang-format off
    options.add_options()
    ("nHmcSteps", po::value<int>(&hmcsteps)->default_value(10),"The number of HMC steps (i.e. the number of configuration updates in the Markov chain).")
    // this is the optimal value...
    ("useGaugeOnly", po::value<bool>(&use_gauge_only)->default_value(false),"Whether to simulate pure gauge theory with HMC. In this case 'nTimeScales' has to be set to 1.")
    ("useMP", po::value<bool>(&use_mp)->default_value(false),"Whether to use the Mass Preconditioning trick.");
    // clang-format on
}
