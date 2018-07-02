/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2018 Alessandro Sciarra
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

#include "parametersHeatbath.hpp"

int meta::ParametersHeatbath::get_thermalizationsteps() const noexcept
{
    return thermalizationsteps;
}
int meta::ParametersHeatbath::get_heatbathsteps() const noexcept
{
    return heatbathsteps;
}
int meta::ParametersHeatbath::get_overrelaxsteps() const noexcept
{
    return overrelaxsteps;
}
int meta::ParametersHeatbath::get_xi() const noexcept
{
    return xi;
}

bool meta::ParametersHeatbath::get_use_aniso() const noexcept
{
    return use_aniso;
}

meta::ParametersHeatbath::ParametersHeatbath() : options("Heatbath options")
{
    // clang-format off
    options.add_options()
    //todo: this is also needed in the HMC!
    ("nThermalizationSteps", po::value<int>(&thermalizationsteps)->default_value(0), "The number of thermalization steps (for the SU(3) Heatbath executable).")
    ("nHeatbathSteps", po::value<int>(&heatbathsteps)->default_value(1000),"The number of heat bath steps (i.e. the number of configuration updates in the Markov chain).")
    ("nOverrelaxationSteps", po::value<int>(&overrelaxsteps)->default_value(1),"The number of overrelaxation steps in the update of a gaugefield configuration.")
    ("useAnisotropy", po::value<bool>(&use_aniso)->default_value(false), "Whether to use an anisotropic lattice, having a lattice spacing different in time and in space directions.")
    ("xi", po::value<int>(&xi)->default_value(1), "The anisotropy coefficient.");
    // clang-format on
}
