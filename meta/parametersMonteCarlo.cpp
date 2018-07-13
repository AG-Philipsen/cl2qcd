/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
 * Copyright (c) 2017,2017,2018 Alessandro Sciarra
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

#include "parametersMonteCarlo.hpp"

int meta::ParametersMonteCarlo::get_thermalizationsteps() const noexcept
{
    return thermalizationsteps;
}
int meta::ParametersMonteCarlo::get_heatbathsteps() const noexcept
{
    return heatbathsteps;
}
int meta::ParametersMonteCarlo::get_overrelaxsteps() const noexcept
{
    return overrelaxsteps;
}
int meta::ParametersMonteCarlo::get_xi() const noexcept
{
    return xi;
}
bool meta::ParametersMonteCarlo::get_use_aniso() const noexcept
{
    return use_aniso;
}
int meta::ParametersMonteCarlo::get_hmcsteps() const noexcept
{
    return hmcsteps;
}
bool meta::ParametersMonteCarlo::get_use_gauge_only() const noexcept
{
    return use_gauge_only;
}

bool meta::ParametersMonteCarlo::get_use_mp() const noexcept
{
    return use_mp;
}
double meta::ParametersMonteCarlo::get_num_tastes() const noexcept
{
    return numberOfTastes;
}
unsigned int meta::ParametersMonteCarlo::get_num_tastes_decimal_digits() const noexcept
{
    return numberOfDecimalDigitsInNumberOfTastes;
}
unsigned int meta::ParametersMonteCarlo::get_num_pseudofermions() const noexcept
{
    return numberOfPseudofermions;
}
unsigned int meta::ParametersMonteCarlo::get_rhmcsteps() const noexcept
{
    return numberOfRhmcSteps;
}

meta::ParametersMonteCarlo::ParametersMonteCarlo()
    : thermalizationsteps(0)
    , heatbathsteps(1000)
    , hmcsteps(10)
    , numberOfRhmcSteps(10)
    , overrelaxsteps(1)
    , xi(1)
    , use_aniso(false)
    , use_gauge_only(false)
    , use_mp(false)
    , numberOfTastes(2)
    , numberOfDecimalDigitsInNumberOfTastes(0)
    , numberOfPseudofermions(1)
    , options("Monte Carlo options")
{
    // clang-format off
    options.add_options()
    ("nThermalizationSteps", po::value<int>(&thermalizationsteps)->default_value(thermalizationsteps), "The number of thermalization steps (for the SU(3) Heatbath executable).")
    ("nHeatbathSteps", po::value<int>(&heatbathsteps)->default_value(heatbathsteps),"The number of heat bath steps (i.e. the number of configuration updates in the Markov chain).")
    ("nOverrelaxationSteps", po::value<int>(&overrelaxsteps)->default_value(overrelaxsteps),"The number of overrelaxation steps in the update of a gaugefield configuration.")
    ("useAnisotropy", po::value<bool>(&use_aniso)->default_value(use_aniso), "Whether to use an anisotropic lattice, having a lattice spacing different in time and in space directions.")
    ("xi", po::value<int>(&xi)->default_value(xi), "The anisotropy coefficient.")
    ("nHmcSteps", po::value<int>(&hmcsteps)->default_value(hmcsteps),"The number of HMC steps (i.e. the number of configuration updates in the Markov chain).")
    ("useGaugeOnly", po::value<bool>(&use_gauge_only)->default_value(use_gauge_only),"Whether to simulate pure gauge theory with HMC. In this case 'nTimeScales' has to be set to 1.")
    ("useMP", po::value<bool>(&use_mp)->default_value(use_mp),"Whether to use the Mass Preconditioning trick.")
    ("nTastes", po::value<double>(&numberOfTastes)->default_value(numberOfTastes, meta::getDefaultForHelper(numberOfTastes)), "The number of tastes of staggered fermions.")
    ("nTastesDecimalDigits", po::value<unsigned int>(&numberOfDecimalDigitsInNumberOfTastes)->default_value(numberOfDecimalDigitsInNumberOfTastes), "The number of decimal digits in the number of staggered tastes.")
    ("nPseudoFermions", po::value<unsigned int>(&numberOfPseudofermions)->default_value(numberOfPseudofermions), "The number of pseudo-fermion species in the multiple pseudofermion technique for staggered fermions only.")
    ("nRhmcSteps", po::value<unsigned int>(&numberOfRhmcSteps)->default_value(numberOfRhmcSteps), "The number of RHMC steps (i.e. the number of configuration updates in the Markov chain).");
    // clang-format on
}
