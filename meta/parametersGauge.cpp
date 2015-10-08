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

#include "parametersGauge.hpp"

//gaugefield parameters
double meta::ParametersGauge::get_beta() const noexcept
{
	return beta;
}
double meta::ParametersGauge::get_rho() const noexcept
{
	return rho;
}
int meta::ParametersGauge::get_rho_iter() const noexcept
{
	return rho_iter;
}
common::action meta::ParametersGauge::get_gaugeact() const noexcept
{
	return gaugeact;
}

bool meta::ParametersGauge::get_use_smearing() const noexcept
{
	return use_smearing;
}

meta::ParametersGauge::ParametersGauge()
	: options("Gaugefield options")
{
	options.add_options()
	("beta", po::value<double>(&beta)->default_value(4.0))
	("use_smearing", po::value<bool>(&use_smearing)->default_value(false))
	("rho", po::value<double>(&rho)->default_value(0.))
	("rho_iter", po::value<int>(&rho_iter)->default_value(0))
	("gaugeact", po::value<std::string>()->default_value("wilson"));
}

meta::ParametersGauge::~ParametersGauge() = default;

po::options_description & meta::ParametersGauge::getOptions()
{
	return options;
}
