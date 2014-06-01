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

po::options_description meta::ParametersRhmc::getOptions()
{
	po::options_description options("RHMC options");
	options.add_options()
		("md_approx_ord", po::value<int>(&md_approx_ord)->default_value(8))
		("metro_approx_ord", po::value<int>(&metro_approx_ord)->default_value(15))
		("findminmax_max", po::value<int>(&findminmax_max)->default_value(5000))
		("findminmax_iteration_block_size", po::value<int>(&findminmax_iteration_block_size)->default_value(25), "find_minmax will check the residual only every N iterations")
		("findminmax_prec", po::value<double>(&findminmax_prec)->default_value(1.e-3))
		("conservative", po::value<bool>(&conservative)->default_value(false))
		("num_tastes", po::value<int>(&num_tastes)->default_value(2))
		("approx_lower", po::value<double>(&approx_lower)->default_value(1.e-5))
		("approx_upper", po::value<double>(&approx_upper)->default_value(1.))
		("rhmcsteps", po::value<int>(&rhmcsteps)->default_value(10))
		("approx_heatbath_file", po::value<std::string>(&approx_heatbath_file)->default_value("Approx_Heatbath"))
		("approx_md_file", po::value<std::string>(&approx_md_file)->default_value("Approx_MD"))
		("approx_metropolis_file", po::value<std::string>(&approx_metropolis_file)->default_value("Approx_Metropolis"))
		("read_rational_approximations_from_file", po::value<bool>(&read_rational_approximations_from_file)->default_value(true));
	return options;
}
