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

#include "parametersHmc.hpp"

#include "../host_functionality/logger.hpp"
#include <stdexcept>

using namespace meta;

double ParametersHmc::get_tau() const noexcept
{
	return tau;
}
bool ParametersHmc::get_reversibility_check() const noexcept
{
	return reversibility_check;
}
int ParametersHmc::get_integrationsteps(size_t timescale) const noexcept
{
	switch(timescale) 
	{
	case 0:
		return integrationsteps0;
	case 1:
		return integrationsteps1;
	case 2:
		return integrationsteps2;
	default:
		throw std::out_of_range("No such timescale");
	}
}
int ParametersHmc::get_hmcsteps() const noexcept
{
	return hmcsteps;
}
int ParametersHmc::get_num_timescales() const noexcept
{
	return num_timescales;
}
ParametersHmc::integrator ParametersHmc::get_integrator(size_t timescale) const noexcept
{
	switch(timescale) 
	{
	case 0:
		return integrator0;
	case 1:
		return integrator1;
	case 2:
		return integrator2;
	default:
		throw std::out_of_range("No such timescale");
	}
}
double ParametersHmc::get_lambda(size_t timescale) const noexcept
{
	switch(timescale)
	{
	case 0:
		return lambda0;
	case 1:
		return lambda1;
	case 2:
		return lambda2;
	default:
		throw std::out_of_range("No such timescale");
	}
}
