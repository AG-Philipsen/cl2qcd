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

#include "parametersObs.hpp"

bool meta::ParametersObs::get_measure_transportcoefficient_kappa() const noexcept
{
	return measure_transportcoefficient_kappa;
}

bool meta::ParametersObs::get_measure_rectangles() const noexcept
{
	return measure_rectangles;
}

bool meta::ParametersObs::get_measure_correlators() const noexcept
{
	return measure_correlators;
}

bool meta::ParametersObs::get_measure_pbp() const noexcept
{
	return measure_pbp;
}

common::pbp_version meta::ParametersObs::get_pbp_version() const noexcept
{
	return pbp_version_;
}

int meta::ParametersObs::get_corr_dir() const noexcept
{
	return corr_dir;
}

int meta::ParametersObs::get_pbp_measurements() const noexcept
{
	return pbp_measurements;
}

meta::ParametersObs::ParametersObs()
	: options("Observables options")
{
	options.add_options()
	("corr_dir", po::value<int>(&corr_dir)->default_value(3), "Direction for the correlator")
	("measure_correlators", po::value<bool>(&measure_correlators)->default_value(true), "Measure fermionic correlators")
	("measure_pbp", po::value<bool>(&measure_pbp)->default_value(false), "Measure chiral condensate")
	("pbp_version",  po::value<std::string>()->default_value("std"), "Version of chiral condensate")
	("pbp_measurements", po::value<int>(&pbp_measurements)->default_value(1), "Number of chiral condensate measurements (stagg. only!)")
	("measure_transportcoefficient_kappa", po::value<bool>(&measure_transportcoefficient_kappa)->default_value(false) )
	("measure_rectangles", po::value<bool>(&measure_rectangles)->default_value(false) );
}

meta::ParametersObs::~ParametersObs() = default;

po::options_description & meta::ParametersObs::getOptions()
{
	return options;
}
