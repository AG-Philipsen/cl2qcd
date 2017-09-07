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

#ifndef _META_PARAMETERS_CONFIG_HPP_
#define _META_PARAMETERS_CONFIG_HPP_

#include "parametersBasic.hpp"

namespace meta {
class ParametersConfig {
public:
	size_t get_precision() const noexcept;

	const std::vector<int> get_selected_devices() const noexcept;
	int get_device_count() const noexcept;
	bool get_use_gpu() const noexcept;
	bool get_use_cpu() const noexcept;
	bool get_enable_profiling() const noexcept;
	int get_nspace() const noexcept;
	int get_ntime() const noexcept;

	//should this go into IO?
	common::startcondition get_startcondition() const noexcept;
	std::string get_sourcefile() const noexcept;
	bool get_ignore_checksum_errors() const noexcept;
	bool get_print_to_screen() const noexcept;
	//parameters to read in gauge configurations
	bool get_read_multiple_configs() const noexcept;
	int get_config_read_start() const noexcept;
	int get_config_read_end() const noexcept;
	int get_config_read_incr() const noexcept;

	std::string get_log_level() const noexcept;
	bool get_use_rec12() const noexcept;
	uint32_t get_host_seed() const noexcept;
	std::string get_initial_prng_state() const noexcept;

	bool get_use_same_rnd_numbers() const noexcept;
	bool is_ocl_compiler_opt_disabled() const noexcept;
	bool get_split_cpu() const noexcept;
	int get_benchmarksteps() const noexcept;

	bool get_use_aniso() const noexcept;

private:
	po::options_description options;

	size_t precision;

	std::vector<int> selected_devices;
	int device_count;
	bool use_gpu;
	bool use_cpu;
	bool enable_profiling;

	bool use_aniso;
	int nspace;
	int ntime;

	//parameters to read in gauge configurations
	bool read_multiple_configs;
	int config_read_start;
	int config_read_end;
	int config_read_incr;
	std::string sourcefile;
	bool ignore_checksum_errors;
	bool print_to_screen;
	uint32_t host_seed;
	std::string initial_prng_state;

	int benchmarksteps;
	bool use_same_rnd_numbers;
	bool use_rec12;
	bool ocl_compiler_opt_disabled;

	std::string log_level;
	bool split_cpu;

protected:
	ParametersConfig();
	virtual ~ParametersConfig();
	ParametersConfig(ParametersConfig const&) = delete;
	ParametersConfig & operator=(ParametersConfig const&) = delete;
	po::options_description & getOptions();

	common::startcondition _startcondition;
};

}

#endif
