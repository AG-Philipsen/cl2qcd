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

int ParametersRhmc::get_md_approx_ord() const noexcept
{
	return md_approx_ord;
}
int ParametersRhmc::get_metro_approx_ord() const noexcept
{
	return metro_approx_ord;
}
int ParametersRhmc::get_findminmax_iteration_block_size() const noexcept
{
	return findminmax_iteration_block_size;
}
int ParametersRhmc::get_findminmax_max() const noexcept
{
	return findminmax_max;
}
double ParametersRhmc::get_findminmax_prec() const noexcept
{
	return findminmax_prec;
}
bool ParametersRhmc::get_conservative() const noexcept
{
	return conservative;
}
int ParametersRhmc::get_num_tastes() const noexcept
{
	return num_tastes;
}
double ParametersRhmc::get_approx_lower() const noexcept
{
	return approx_lower;
}
double ParametersRhmc::get_approx_upper() const noexcept
{
	return approx_upper;
}
int ParametersRhmc::get_rhmcsteps() const noexcept
{
	return rhmcsteps;
}
std::string ParametersRhmc::get_approx_heatbath_file() const noexcept
{
	return approx_heatbath_file;
}
std::string ParametersRhmc::get_approx_md_file() const noexcept
{
	return approx_md_file;
}
std::string ParametersRhmc::get_approx_metropolis_file() const noexcept
{
	return approx_metropolis_file;
}
bool ParametersRhmc::get_read_rational_approximations_from_file() const noexcept
{
	return read_rational_approximations_from_file;
}
