/** @file
 * Generic utility functions
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#ifndef _META_UTIL_
#define _META_UTIL_

#include "inputparameters.hpp"

#include "../common_header_files/globaldefs.h"
#include "../common_header_files/types.h"
#include <iostream>
#include <sstream>
#include <string.h>
#include <utility>

namespace meta {
size_t get_volspace(const Inputparameters&);
size_t get_volspace(const int ns);
size_t get_vol4d(const Inputparameters&);
size_t get_vol4d(const int nt, const int ns);
bool get_use_rectangles(const Inputparameters& params);
hmc_float get_mubar(const Inputparameters& params);
hmc_float get_mubar_mp(const Inputparameters& params);
size_t get_float_size(const Inputparameters& params);
size_t get_mat_size(const Inputparameters& params);
size_t get_plaq_norm(const Inputparameters& params);
size_t get_tplaq_norm(const Inputparameters& params);
size_t get_splaq_norm(const Inputparameters& params);
size_t get_rect_norm(const Inputparameters& params);
size_t get_poly_norm(const Inputparameters& params);
double get_c0(const Inputparameters& params);
double get_c1(const Inputparameters& params);
double get_xi_0(const Inputparameters& params);

void print_info_global(const meta::Inputparameters& params);
void print_info_global(std::ostream* os, const meta::Inputparameters& params);
void print_info_configs_io(const meta::Inputparameters& params);
void print_info_configs_io(std::ostream * os, const meta::Inputparameters& params);
void print_info_prng_io(const meta::Inputparameters& params);
void print_info_prng_io(std::ostream * os, const meta::Inputparameters& params);
void print_info_observables_gauge_io(const meta::Inputparameters& params);
void print_info_observables_gauge_io(std::ostream * os, const meta::Inputparameters& params);
void print_info_heatbath(const Inputparameters& params);
void print_info_heatbath(std::ostream* os, const Inputparameters& params);
void print_info_inverter(const Inputparameters& params);
void print_info_inverter(std::ostream* os, const Inputparameters& params);
void print_info_hmc(const Inputparameters& params);
void print_info_hmc(std::ostream* os, const Inputparameters& params);
void print_info_rhmc(const Inputparameters& params);
void print_info_rhmc(std::ostream* os, const Inputparameters& params);
std::string get_ferm_obs_corr_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
std::string get_ferm_obs_pbp_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
std::string get_gauge_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
std::string get_hmc_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
std::string get_rhmc_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
std::string create_prng_name(const Inputparameters& parameters, int number) noexcept;
std::string create_profiling_data_filename(const Inputparameters& parameters, std::string executableName) noexcept;
std::string create_prng_name(const Inputparameters& parameters) noexcept;
void print_info_flavour_doublet_correlators(const meta::Inputparameters& params);
void print_info_flavour_doublet_correlators(std::ostream * os, const meta::Inputparameters& params);
template<typename Container> void free_container(Container& c) noexcept;
std::string createLogfileName(const char* name);
std::pair<int, std::vector<const char*>> addOptionsToArgv(int argc, const char** argv, std::vector<const char*> op);

/*
* TEMPLATE IMPLEMENTATIONS
*/
template<typename Container> void free_container(Container& c) noexcept {
	while(!c.empty()) {
		delete c.back();
		c.pop_back();
	}
}


}

#endif /* META_UTIL_ */
