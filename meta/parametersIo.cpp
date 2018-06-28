/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Alessandro Sciarra
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

#include "parametersIo.hpp"

int meta::ParametersIo::get_writefrequency() const noexcept
{
    return writefrequency;
}
int meta::ParametersIo::get_savefrequency() const noexcept
{
    return savefrequency;
}
int meta::ParametersIo::get_savepointfrequency() const noexcept
{
    return savepointfrequency;
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

meta::ParametersIo::ParametersIo() : options("IO options")
{
    // clang-format off
	options.add_options()
	("onlineMeasureEvery", po::value<int>(&writefrequency)->default_value(1), "The frequency of online measurements.")
	("createConfAndPrngCheckpointsEvery", po::value<int>(&savefrequency)->default_value(100), "The frequency at which gaugefield configuration and prng state ckeckpoints (conf.#####, prng.#####) are saved.")
    ("saveConfAndPrngEvery", po::value<int>(&savepointfrequency)->default_value(1), "The frequency at which gaugefield configuration and prng state (conf.save, prng.save) are saved.")
	("nDigitsInConfCheckpoint", po::value<int>(&config_number_digits)->default_value(5), "The number of digits to name gaugefield configuration checkpoints.")
	("confPrefix", po::value<std::string>(&config_prefix)->default_value("conf."), "The prefix for gaugefield configuration file name.")
	("confPostfix", po::value<std::string>(&config_postfix)->default_value(""), "The postfix for gaugefield configuration file name.")
	("prngPrefix", po::value<std::string>(&prng_prefix)->default_value("prng."), "The prefix for prng state file name.")
	("prngPostfix", po::value<std::string>(&prng_postfix)->default_value(""), "The postfix for prng state file name.")
	("rectanglesFilename", po::value<std::string>(&rectanglesFilename)->default_value("gaugeObsRectangles.dat"), "The filename for rectangles measurements.")
	("transportCoefficientKappaFileName", po::value<std::string>(&transportcoefficientKappaFilename)->default_value("GaugeObsKappa"), "The filename for transport coefficient kappa measurements.")
	("profilingDataPrefix", po::value<std::string>(&profiling_data_prefix)->default_value(""), "The prefix for profiling data file name.")
	("profilingDataPostfix", po::value<std::string>(&profiling_data_postfix)->default_value("_profiling_data"), "The postfix for profiling data file name.")
	("gaugeObsInSingleFile", po::value<bool>(&gauge_obs_to_single_file)->default_value(true), "Whether to save gauge observables in a single file. These are plaquette and polyakov measurements performed at the beginning of the 'thermalization' phase even when 'nThermalizationSteps' is set to 0.")
	("gaugeObsPrefix", po::value<std::string>(&gauge_obs_prefix)->default_value("gaugeObs"), "The prefix for the gauge observables file name.")
	("gaugeObsPostfix", po::value<std::string>(&gauge_obs_postfix)->default_value(".dat"), "The postfix for the gauge observables file name.")
	("fermObsInSingleFile", po::value<bool>(&ferm_obs_to_single_file)->default_value(false), "Whether to save fermionic observables (chiral condensate, correlators) in a single file.")
	("fermObsCorrelatorsPrefix", po::value<std::string>(&ferm_obs_corr_prefix)->default_value(""), "The prefix for fermionic observables file name for correlator measurements.")
	("fermObsCorrelatorsPostfix", po::value<std::string>(&ferm_obs_corr_postfix)->default_value("_correlators.dat"), "The postfix for fermionic observables file name for correlator measurements.")
	("fermObsPbpPrefix", po::value<std::string>(&ferm_obs_pbp_prefix)->default_value(""), "The prefix for fermionic observables file name for chiral condensate measurements.")
	("fermObsPbpPostfix", po::value<std::string>(&ferm_obs_pbp_postfix)->default_value("_pbp.dat"), "The postfix for fermionic observables file name for chiral condensate measurements.")
	("hmcObsToSingleFile", po::value<bool>(&hmc_obs_to_single_file)->default_value(true), "Whether to save hmc observables to a single file. These are plaq, tplaq, splaq, poly.re, poly.im, |poly|, deltaH, acceptance, timeTrajectory.")
	("hmcObsPrefix", po::value<std::string>(&hmc_obs_prefix)->default_value("hmc_output"), "The prefix for hmc observables file name.")
	("hmcObsPostfix", po::value<std::string>(&hmc_obs_postfix)->default_value(""), "The postfix for hmc observables file name.")
	("rhmcObsToSingleFile", po::value<bool>(&rhmc_obs_to_single_file)->default_value(true), "Whether to save rhmc observables to one single file. These are the same as for hmc.")
	("rhmcObsPrefix", po::value<std::string>(&rhmc_obs_prefix)->default_value("rhmc_output"), "The prefix for rhmc observables file name.")
	("rhmcObsPostfix", po::value<std::string>(&rhmc_obs_postfix)->default_value(""), "The postfix for rhmc observables file name.");
    // clang-format on
}
