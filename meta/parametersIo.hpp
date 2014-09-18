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

#ifndef _META_PARAMETERS_IO_HPP_
#define _META_PARAMETERS_IO_HPP_

#include "parametersBasic.hpp"

namespace meta {
class ParametersIo {
public:
	int get_writefrequency() const noexcept;
	int get_savefrequency() const noexcept;
	int get_config_number_digits() const noexcept;
	std::string get_prng_prefix() const noexcept;
	std::string get_prng_postfix() const noexcept;
	std::string get_config_prefix() const noexcept;
	std::string get_config_postfix() const noexcept;
	std::string get_ferm_obs_corr_prefix() const noexcept;
	std::string get_ferm_obs_corr_postfix() const noexcept;
	std::string get_ferm_obs_pbp_prefix() const noexcept;
	std::string get_ferm_obs_pbp_postfix() const noexcept;
	std::string get_gauge_obs_prefix() const noexcept;
	std::string get_gauge_obs_postfix() const noexcept;
	bool get_ferm_obs_to_single_file() const noexcept;
	bool get_gauge_obs_to_single_file() const noexcept;
	std::string get_hmc_obs_prefix() const noexcept;
	std::string get_hmc_obs_postfix() const noexcept;
	bool get_hmc_obs_to_single_file() const noexcept;
	std::string get_rhmc_obs_prefix() const noexcept;
	std::string get_rhmc_obs_postfix() const noexcept;
	bool get_rhmc_obs_to_single_file() const noexcept;
	std::string get_profiling_data_prefix() const noexcept;
	std::string get_profiling_data_postfix() const noexcept;
	std::string get_rectanglesFilename() const noexcept;
	std::string get_transportcoefficientKappaFilename() const noexcept;

private:
	po::options_description options;
	int writefrequency;
	int savefrequency;
	int config_number_digits;
	std::string config_prefix;
	std::string config_postfix;
	std::string prng_prefix;
	std::string prng_postfix;
	std::string rectanglesFilename;
	std::string transportcoefficientKappaFilename;
	std::string profiling_data_prefix;
	std::string profiling_data_postfix;
	bool gauge_obs_to_single_file;
	std::string gauge_obs_prefix;
	std::string gauge_obs_postfix;
	bool ferm_obs_to_single_file;
	std::string ferm_obs_corr_prefix;
	std::string ferm_obs_corr_postfix;
	std::string ferm_obs_pbp_prefix;
	std::string ferm_obs_pbp_postfix;
	bool hmc_obs_to_single_file;
	std::string hmc_obs_prefix;
	std::string hmc_obs_postfix;
	bool rhmc_obs_to_single_file;
	std::string rhmc_obs_prefix;
	std::string rhmc_obs_postfix;

protected:
	ParametersIo();
	po::options_description & getOptions();
};

}

#endif
