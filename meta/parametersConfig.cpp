/** @file
 *
 * Copyright (c) 2014,2015 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2018 Alessandro Sciarra
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

#include "parametersConfig.hpp"

#include "../executables/exceptions.hpp"

#include <boost/algorithm/string.hpp>

static common::startcondition translateStartConditionToEnum(std::string);

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

// parameters to read in gauge configurations
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
    : precision(sizeof(double) * 8)
    , selected_devices(0)
    , device_count(0)
    , use_gpu(true)
    , use_cpu(true)
    , enable_profiling(false)
    , nspace(4)
    , ntime(8)
    , read_multiple_configs(false)
    , config_read_start(0)
    , config_read_end(1)
    , config_read_incr(1)
    , sourcefile("conf.00000")
    , ignore_checksum_errors(false)
    , print_to_screen(false)
    , host_seed(4815)
    , initial_prng_state("")
    , benchmarksteps(500)
    , use_same_rnd_numbers(false)
    , use_rec12(false)
    , ocl_compiler_opt_disabled(false)
    , log_level("ALL")
    , split_cpu(false)
    , options("Configuration options")
    , _startconditionString("cold")
    , _startcondition(common::startcondition::cold_start)
{
    // clang-format off
    options.add_options()
    ("confPrecision", po::value<size_t>(&precision)->default_value(precision), "The precision in bit for the storage of gauge configurations in conf files.")
    ("deviceId", po::value<std::vector<int>>(&selected_devices), "The ID number of a device to use. Can be specified multiple times to use multiple devices.")
    ("nDevices", po::value<int>(&device_count)->default_value(device_count), "The maximum number of devices to use.")
    ("useGPU", po::value<bool>(&use_gpu)->default_value(use_gpu), "Whether to use GPUs.")
    ("useCPU", po::value<bool>(&use_cpu)->default_value(use_cpu), "Whether to use CPUs.")
    ("enableProfiling", po::value<bool>(&enable_profiling)->default_value(enable_profiling), "Whether to profile kernel execution. This option implies slower performance due to synchronization after each kernel call.")
    ("nSpace", po::value<int>(&nspace)->default_value(nspace), "The spatial extent of the lattice.")
    ("nTime", po::value<int>(&ntime)->default_value(ntime), "The temporal extent of the lattice.")
    ("startCondition", po::value<std::string>(&_startconditionString)->default_value(_startconditionString), "The gaugefield starting condition (e.g. cold, hot, continue).")
    ("initialConf", po::value<std::string>(&sourcefile)->default_value(sourcefile), "The path of the file containing the gauge configuration to start from.")
    ("printToScreen", po::value<bool>(&print_to_screen)->default_value(print_to_screen), "Whether to print the onfly measurements to the standard output.")
    ("hostSeed", po::value<uint32_t>(&host_seed)->default_value(host_seed),"The random seed to initialize the pseudo random number generator.")
    ("initialPRNG", po::value<std::string>(&initial_prng_state)->default_value(initial_prng_state),"The path of the file containing the pseudo random number generator state to start from.")
    ("useSameRandomNumbers", po::value<bool>(&use_same_rnd_numbers)->default_value(use_same_rnd_numbers), "Whether to use random numbers compatible with a scalar version. If it is set to 'true', then the number of random states is one instead of being equal to the number of global threads. This option implies a huge loss in performance, hence it should be used with care!")
    ("disableOclCompilerOptimization", po::value<bool>(&ocl_compiler_opt_disabled)->default_value(ocl_compiler_opt_disabled), "Whether to disable OpenCL compiler from performing optimizations (cf. -cl-disable-opt).")
    ("useReconstruct12", po::value<bool>(&use_rec12)->default_value(use_rec12), "Whether to use reconstruct 12 compression for SU3 matrices, i.e. consider gauge links stored using 12 complex numbers instead of 18.")
    ("logLevel", po::value<std::string>(&log_level)->default_value(log_level), "The minimum output log level (one among ALL, TRACE, DEBUG, INFO, WARN, ERROR, FATAL, OFF).")
    ("readMultipleConfs", po::value<bool>(&read_multiple_configs)->default_value(read_multiple_configs), "Whether to use more than one gaugefield configuration at once.")
    ("readFromConfNumber", po::value<int>(&config_read_start)->default_value(config_read_start), "The number to begin with when using more than one gaugefield configuration at once.")
    ("readUntilConfNumber", po::value<int>(&config_read_end)->default_value(config_read_end), "The number to end at when using more than one gaugefield configuration at once.")
    ("readConfsEvery", po::value<int>(&config_read_incr)->default_value(config_read_incr), "The increment for the gaugefield configuration number when using more than one gaugefield configuration at once.")
    ("splitCPU", po::value<bool>(&split_cpu)->default_value(split_cpu), "Whether to split the CPU into multiple devices to avoid numa issues. This option requires OpenCL 1.2 at least.")
    ("benchmarkIterations", po::value<int>(&benchmarksteps)->default_value(benchmarksteps), "The number of times a kernel is executed for benchmark purposes.")
    ("ignoreChecksumErrors", po::value<bool>(&ignore_checksum_errors)->default_value(ignore_checksum_errors), "Whether to ignore checksum errors, e.g. reading conf files.");
    // clang-format on
}

static common::startcondition translateStartConditionToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::startcondition> m;
    m["cold_start"]        = common::cold_start;
    m["cold"]              = common::cold_start;
    m["hot_start"]         = common::hot_start;
    m["hot"]               = common::hot_start;
    m["start_from_source"] = common::start_from_source;
    m["source"]            = common::start_from_source;
    m["continue"]          = common::start_from_source;

    common::startcondition a = m[s];
    if (a) {
        return a;
    } else {
        throw Invalid_Parameters("Unkown start condition!",
                                 "cold_start, cold, hot_start, hot, start_from_source, source, continue", s);
    }
}

void meta::ParametersConfig::makeNeededTranslations()
{
    _startcondition = translateStartConditionToEnum(_startconditionString);
}
