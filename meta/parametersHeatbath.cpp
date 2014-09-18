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

meta::ParametersHeatbath::ParametersHeatbath()
	: options("Heatbath options")
{
	options.add_options()
	//todo: this is also needed in the HMC!
	("thermalization", po::value<int>(&thermalizationsteps)->default_value(0))
	("heatbathsteps", po::value<int>(&heatbathsteps)->default_value(1000))
	("overrelaxsteps", po::value<int>(&overrelaxsteps)->default_value(1))
	//todo: is this used?
	("xi", po::value<int>(&xi)->default_value(1));
}

po::options_description & meta::ParametersHeatbath::getOptions()
{
	return options;
}
