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

#include "parametersConfig.hpp"

size_t meta::ParametersConfig::get_precision() const noexcept
{
	return precision;
}

const std::vector<int> meta::ParametersConfig::get_selected_devices() const noexcept
{
	return selected_devices;
}
int meta::ParametersConfig::get_device_count() const noexcept
{
	return device_count;
}
bool meta::ParametersConfig::get_use_gpu() const noexcept
{
	return use_gpu;
}
bool meta::ParametersConfig::get_use_cpu() const noexcept
{
	return use_cpu;
}
bool meta::ParametersConfig::get_enable_profiling() const noexcept
{
	return enable_profiling;
}

bool meta::ParametersConfig::get_use_aniso() const noexcept
{
	return use_aniso;
}

int meta::ParametersConfig::get_nspace() const noexcept
{
	return nspace;
}
int meta::ParametersConfig::get_ntime() const noexcept
{
	return ntime;
}

common::startcondition meta::ParametersConfig::get_startcondition() const noexcept
{
	return _startcondition;
}

std::string meta::ParametersConfig::get_sourcefile() const noexcept
{
	return sourcefile;
}
bool meta::ParametersConfig::get_ignore_checksum_errors() const noexcept
{
	return ignore_checksum_errors;
}
bool meta::ParametersConfig::get_print_to_screen() const noexcept
{
	return print_to_screen;
}
uint32_t meta::ParametersConfig::get_host_seed() const noexcept
{
	return host_seed;
}
std::string meta::ParametersConfig::get_initial_prng_state() const noexcept
{
	return initial_prng_state;
}

int meta::ParametersConfig::get_benchmarksteps() const noexcept
{
	return benchmarksteps;
}

bool meta::ParametersConfig::get_use_same_rnd_numbers() const noexcept
{
	return use_same_rnd_numbers;
}

bool meta::ParametersConfig::is_ocl_compiler_opt_disabled() const noexcept
{
	return ocl_compiler_opt_disabled;
}

std::string meta::ParametersConfig::get_log_level() const noexcept
{
	return log_level;
}

//parameters to read in gauge configurations
bool meta::ParametersConfig::get_read_multiple_configs() const noexcept
{
	return read_multiple_configs;
}
int meta::ParametersConfig::get_config_read_start() const noexcept
{
	return config_read_start;
}
int meta::ParametersConfig::get_config_read_end() const noexcept
{
	return config_read_end;
}
int meta::ParametersConfig::get_config_read_incr() const noexcept
{
	return config_read_incr;
}

bool meta::ParametersConfig::get_use_rec12() const noexcept
{
	return use_rec12;
}

bool meta::ParametersConfig::get_split_cpu() const noexcept
{
	return split_cpu;
}

meta::ParametersConfig::ParametersConfig()
	: options("Configuration options")
{
	options.add_options()
	("prec", po::value<size_t>(&precision)->default_value(sizeof(double) * 8))
	("device,d", po::value<std::vector<int>>(&selected_devices), "ID of a device to use. Can be specified multiple times.")
	("num_dev", po::value<int>(&device_count)->default_value(0), "Maximum number of devices to use.")
	("use_gpu", po::value<bool>(&use_gpu)->default_value(true), "Use GPUs")
	("use_cpu", po::value<bool>(&use_cpu)->default_value(true), "Use CPUs")
	("enable_profiling", po::value<bool>(&enable_profiling)->default_value(false), "Enable profiling of kernel execution. Implies slower performance due to synchronization after each kernel call.")
	("nspace", po::value<int>(&nspace)->default_value(4))
	("ntime", po::value<int>(&ntime)->default_value(8))
	("startcondition", po::value<std::string>()->default_value("cold_start"))
	("sourcefile", po::value<std::string>(&sourcefile)->default_value("conf.00000"))
	("print_to_screen", po::value<bool>(&print_to_screen)->default_value(false))
	("host_seed", po::value<uint32_t>(&host_seed)->default_value(4815))
	("initial_prng_state", po::value<std::string>(&initial_prng_state)->default_value(""))
	("use_same_rnd_numbers", po::value<bool>(&use_same_rnd_numbers)->default_value(false), "Use random numbers compatible with a scalar version. SLOW!")
	("disable-ocl-compiler-opt", po::value<bool>(&ocl_compiler_opt_disabled)->default_value(false), "Disable OpenCL compiler from performing optimizations (adds -cl-disable-opt)")
	("use_rec12", po::value<bool>(&use_rec12)->default_value(false), "Use reconstruct 12 compression for SU3 matrices")
	("log-level", po::value<std::string>(&log_level)->default_value("ALL"), "Minimum output log level: ALL TRACE DEBUG INFO WARN ERROR FATAL OFF")
	("read_multiple_configs", po::value<bool>(&read_multiple_configs)->default_value(false), "Read in more than one gaugefield configuration")
	("config_read_start", po::value<int>(&config_read_start)->default_value(0), "Number to begin with when reading in more than one gaugefield configuration")
	("config_read_end", po::value<int>(&config_read_end)->default_value(1), "Number to end with when reading in more than one gaugefield configuration")
	("config_read_incr", po::value<int>(&config_read_incr)->default_value(1), "Increment for gaugefield configuration number when reading in more than one gaugefield configuration")
	("split_cpu", po::value<bool>(&split_cpu)->default_value(false), "Split the CPU into multiple devices to avoid numa issues. (Requires OpenCL 1.2 at least)")
	("benchmarksteps", po::value<int>(&benchmarksteps)->default_value(500))
	("ignore_checksum_errors", po::value<bool>(&ignore_checksum_errors)->default_value(false))
	//todo: this is not used ?!
	("use_aniso", po::value<bool>(&use_aniso)->default_value(false));
}

meta::ParametersConfig::~ParametersConfig() = default;

po::options_description & meta::ParametersConfig::getOptions()
{
	return options;
}
