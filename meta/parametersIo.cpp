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

#include "parametersIo.hpp"

int meta::ParametersIo::get_writefrequency() const noexcept
{
	return writefrequency;
}
int meta::ParametersIo::get_savefrequency() const noexcept
{
	return savefrequency;
}

int meta::ParametersIo::get_config_number_digits() const noexcept
{
	return config_number_digits;
}

std::string meta::ParametersIo::get_profiling_data_prefix() const noexcept
{
	return profiling_data_prefix;
}

std::string meta::ParametersIo::get_profiling_data_postfix() const noexcept
{
	return profiling_data_postfix;
}

std::string meta::ParametersIo::get_prng_prefix() const noexcept
{
	return prng_prefix;
}

std::string meta::ParametersIo::get_prng_postfix() const noexcept
{
	return prng_postfix;
}

std::string meta::ParametersIo::get_config_prefix() const noexcept
{
	return config_prefix;
}

std::string meta::ParametersIo::get_config_postfix() const noexcept
{
	return config_postfix;
}

std::string meta::ParametersIo::get_ferm_obs_corr_prefix() const noexcept
{
	return ferm_obs_corr_prefix;
}

std::string meta::ParametersIo::get_ferm_obs_corr_postfix() const noexcept
{
	return ferm_obs_corr_postfix;
}

std::string meta::ParametersIo::get_ferm_obs_pbp_prefix() const noexcept
{
	return ferm_obs_pbp_prefix;
}

std::string meta::ParametersIo::get_ferm_obs_pbp_postfix() const noexcept
{
	return ferm_obs_pbp_postfix;
}

std::string meta::ParametersIo::get_gauge_obs_prefix() const noexcept
{
	return gauge_obs_prefix;
}

std::string meta::ParametersIo::get_gauge_obs_postfix() const noexcept
{
	return gauge_obs_postfix;
}

bool meta::ParametersIo::get_ferm_obs_to_single_file() const noexcept
{
	return ferm_obs_to_single_file;
}

bool meta::ParametersIo::get_gauge_obs_to_single_file() const noexcept
{
	return gauge_obs_to_single_file;
}

std::string meta::ParametersIo::get_hmc_obs_prefix() const noexcept
{
	return hmc_obs_prefix;
}

std::string meta::ParametersIo::get_hmc_obs_postfix() const noexcept
{
	return hmc_obs_postfix;
}

bool meta::ParametersIo::get_hmc_obs_to_single_file() const noexcept
{
	return hmc_obs_to_single_file;
}

std::string meta::ParametersIo::get_rhmc_obs_prefix() const noexcept
{
	return rhmc_obs_prefix;
}

std::string meta::ParametersIo::get_rhmc_obs_postfix() const noexcept
{
	return rhmc_obs_postfix;
}

bool meta::ParametersIo::get_rhmc_obs_to_single_file() const noexcept
{
	return rhmc_obs_to_single_file;
}

std::string meta::ParametersIo::get_rectanglesFilename() const noexcept
{
	return rectanglesFilename;
}
std::string meta::ParametersIo::get_transportcoefficientKappaFilename() const noexcept
{
	return transportcoefficientKappaFilename;
}

po::options_description meta::ParametersIo::getOptions()
{
	po::options_description options("IO options");
	options.add_options()
	("writefrequency", po::value<int>(&writefrequency)->default_value(1))
	("savefrequency", po::value<int>(&savefrequency)->default_value(100))
	("config_number_digits", po::value<int>(&config_number_digits)->default_value(5), "Number of digits to name gaugefield configurations")
	("config_prefix", po::value<std::string>(&config_prefix)->default_value("conf."), "Prefix for gaugefield configuration")
	("config_postfix", po::value<std::string>(&config_postfix)->default_value(""), "Postfix for gaugefield configuration")
	("prng_prefix", po::value<std::string>(&prng_prefix)->default_value("prng."), "Prefix for PRNG configuration")
	("prng_postfix", po::value<std::string>(&prng_postfix)->default_value(""), "Postfix for PRNG configuration")
	("rectanglesFilename", po::value<std::string>(&rectanglesFilename)->default_value("gaugeObsRectangles.dat"), "Filename for rectangles measurements")
	("transportcoefficientKappaFilename", po::value<std::string>(&transportcoefficientKappaFilename)->default_value("GaugeObsKappa"), "Filename for transportcoefficient kappa measurements")
	("profiling_data_prefix", po::value<std::string>(&profiling_data_prefix)->default_value(""), "Prefix for profiling data filename")
	("profiling_data_postfix", po::value<std::string>(&profiling_data_postfix)->default_value("_profiling_data"), "Postfix for profiling data filename")
	("gauge_obs_to_single_file", po::value<bool>(&gauge_obs_to_single_file)->default_value(true), "Save gauge observables to one single file")
	("gauge_obs_prefix", po::value<std::string>(&gauge_obs_prefix)->default_value("gaugeObs"), "Prefix for gauge observables file")
	("gauge_obs_postfix", po::value<std::string>(&gauge_obs_postfix)->default_value(".dat"), "Postfix for gauge observables file")
	("ferm_obs_to_single_file", po::value<bool>(&ferm_obs_to_single_file)->default_value(false), "Save fermionic observables to one single file")
	("ferm_obs_corr_prefix", po::value<std::string>(&ferm_obs_corr_prefix)->default_value(""), "Prefix for fermionic observables (correlators) file")
	("ferm_obs_corr_postfix", po::value<std::string>(&ferm_obs_corr_postfix)->default_value("_correlators.dat"), "Postfix for fermionic observables (correlators) file")
	("ferm_obs_pbp_prefix", po::value<std::string>(&ferm_obs_pbp_prefix)->default_value(""), "Prefix for fermionic observables (chiral condensate) file")
	("ferm_obs_pbp_postfix", po::value<std::string>(&ferm_obs_pbp_postfix)->default_value("_pbp.dat"), "Postfix for fermionic observables (chiral condensate) file")
	("hmc_obs_to_single_file", po::value<bool>(&hmc_obs_to_single_file)->default_value(true), "Save hmc observables to one single file")
	("hmc_obs_prefix", po::value<std::string>(&hmc_obs_prefix)->default_value("hmc_output"), "Prefix for hmc observables file")
	("hmc_obs_postfix", po::value<std::string>(&hmc_obs_postfix)->default_value(""), "Postfix for hmc observables file")
	("rhmc_obs_to_single_file", po::value<bool>(&rhmc_obs_to_single_file)->default_value(true), "Save rhmc observables to one single file")
	("rhmc_obs_prefix", po::value<std::string>(&rhmc_obs_prefix)->default_value("rhmc_output"), "Prefix for rhmc observables file")
	("rhmc_obs_postfix", po::value<std::string>(&rhmc_obs_postfix)->default_value(""), "Postfix for rhmc observables file");

	return options;
}

