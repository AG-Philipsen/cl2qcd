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

#include "util.hpp"

size_t meta::get_volspace(const Inputparameters& params)
{
	return meta::get_volspace(params.get_nspace());
}
size_t meta::get_volspace(const int ns)
{
	return ns * ns * ns;
}

size_t meta::get_vol4d(const Inputparameters& params)
{
	return meta::get_vol4d(params.get_ntime(), params.get_nspace());
}
size_t meta::get_vol4d(const int nt, const int ns)
{
    return meta::get_volspace(ns) * nt;
}

bool meta::get_use_rectangles(const Inputparameters& params)
{
	switch(params.get_gaugeact()) {
		case common::action::tlsym:
		case common::action::iwasaki:
		case common::action::dbw2:
			return true;
		default:
			return false;
	}
}
hmc_float meta::get_mubar(const Inputparameters& params)
{
	return 2. * params.get_kappa() * params.get_mu();
}
hmc_float meta::get_mubar_mp(const Inputparameters& params)
{
	return 2. * params.get_kappa_mp() * params.get_mu_mp();
}

size_t meta::get_float_size(const Inputparameters& params)
{
	return params.get_precision() / 8;
}
size_t meta::get_mat_size(const Inputparameters&)
{
	// TODO with rec12 this becomes 6
	return 9;
}

size_t meta::get_plaq_norm(const Inputparameters& params)
{
	return (get_vol4d(params) * NDIM * (NDIM - 1) ) / 2.;
}
size_t meta::get_tplaq_norm(const Inputparameters& params)
{
	return (get_vol4d(params) * (NDIM - 1));
}
size_t meta::get_splaq_norm(const Inputparameters& params)
{
	return (get_vol4d(params) * (NDIM - 1) * (NDIM - 2)) / 2. ;
}
size_t meta::get_rect_norm(const Inputparameters& params)
{
	return NDIM * (NDIM - 1) * get_vol4d(params);
}
size_t meta::get_poly_norm(const Inputparameters& params)
{
	return get_volspace(params);
}

double meta::get_c0(const Inputparameters& params)
{
	switch(params.get_gaugeact()) {
		case common::action::tlsym:
		case common::action::iwasaki:
		case common::action::dbw2:
			return 1. - 8. * get_c1(params);
		default:
			return 1.;
	}
}
double meta::get_c1(const Inputparameters& params)
{
	switch(params.get_gaugeact()) {
		case common::action::tlsym:
			return -0.083333333;
		case common::action::iwasaki:
			return -0.331;
		case common::action::dbw2:
			return -1.4069;
		default:
			return 0.;
	}
}
double meta::get_xi_0(const Inputparameters& params)
{
	double aniso = params.get_xi();
	double beta = params.get_beta();
	double eta = (1.002503 * aniso * aniso * aniso + .39100 * aniso * aniso + 1.47130 * aniso - 0.19231) /
	             (aniso * aniso * aniso + 0.26287 * aniso * aniso + 1.59008 * aniso - 0.18224);
	return aniso / (1. + (1. - 1. / aniso) * eta / 6. * (1 - 0.55055 * 2 * NC / beta) / (1 - 0.77810 * 2 * NC / beta) * 2 * NC / beta );
}

std::string meta::get_ferm_obs_pbp_file_name(const Inputparameters& parameters, std::string conf_name) noexcept {
	if(parameters.get_ferm_obs_to_single_file() ) {
		return parameters.get_ferm_obs_pbp_prefix() +  parameters.get_ferm_obs_pbp_postfix();
	} else {
		return  parameters.get_ferm_obs_pbp_prefix() + conf_name + parameters.get_ferm_obs_pbp_postfix();
	}
}

std::string meta::get_ferm_obs_corr_file_name(const Inputparameters& parameters, std::string conf_name) noexcept {
	if(parameters.get_ferm_obs_to_single_file() ) {
		return parameters.get_ferm_obs_corr_prefix() +  parameters.get_ferm_obs_corr_postfix();
	} else {
		switch ( parameters.get_startcondition() ) {
			case common::start_from_source :
				return  parameters.get_ferm_obs_corr_prefix() + conf_name + parameters.get_ferm_obs_corr_postfix();
				break;
			case common::hot_start :
				return  parameters.get_ferm_obs_corr_prefix() + "conf.hot" +   parameters.get_ferm_obs_corr_postfix() ;
				break;
			case common::cold_start :
				return  parameters.get_ferm_obs_corr_prefix() + "conf.cold" + parameters.get_ferm_obs_corr_postfix() ;
				break;
			default:
				throw std::invalid_argument("Unknown start condition selected. Don't know how to generate proper observable file name.");
		}
	}
}

std::string meta::get_gauge_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept {
	if(parameters.get_gauge_obs_to_single_file() ) {
		return parameters.get_gauge_obs_prefix() +  parameters.get_gauge_obs_postfix();
	} else {
		switch ( parameters.get_startcondition() ) {
			case common::start_from_source :
				return  parameters.get_gauge_obs_prefix() + conf_name + parameters.get_gauge_obs_postfix();
				break;
			case common::hot_start :
				return  parameters.get_gauge_obs_prefix() + "conf.hot" +   parameters.get_gauge_obs_postfix() ;
				break;
			case common::cold_start :
				return  parameters.get_gauge_obs_prefix() + "conf.cold" + parameters.get_gauge_obs_postfix() ;
				break;
			default:
				throw std::invalid_argument("Unknown start condition selected. Don't know how to generate proper observable file name.");
		}
	}
}

std::string meta::get_hmc_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept {
	if(parameters.get_hmc_obs_to_single_file() ) {
		return parameters.get_hmc_obs_prefix() +  parameters.get_hmc_obs_postfix();
	} else {
		switch ( parameters.get_startcondition() ) {
			case common::start_from_source :
				return  parameters.get_hmc_obs_prefix() + conf_name + parameters.get_hmc_obs_postfix();
				break;
			case common::hot_start :
				return  parameters.get_hmc_obs_prefix() + "conf.hot" +   parameters.get_hmc_obs_postfix() ;
				break;
			case common::cold_start :
				return  parameters.get_hmc_obs_prefix() + "conf.cold" + parameters.get_hmc_obs_postfix() ;
				break;
			default:
				throw std::invalid_argument("Unknown start condition selected. Don't know how to generate proper observable file name.");
		}
	}
}

std::string meta::get_rhmc_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept {
	if(parameters.get_rhmc_obs_to_single_file() ) {
		return parameters.get_rhmc_obs_prefix() +  parameters.get_rhmc_obs_postfix();
	} else {
		switch ( parameters.get_startcondition() ) {
			case common::start_from_source :
				return  parameters.get_rhmc_obs_prefix() + conf_name + parameters.get_rhmc_obs_postfix();
				break;
			case common::hot_start :
				return  parameters.get_rhmc_obs_prefix() + "conf.hot" +   parameters.get_rhmc_obs_postfix() ;
				break;
			case common::cold_start :
				return  parameters.get_rhmc_obs_prefix() + "conf.cold" + parameters.get_rhmc_obs_postfix() ;
				break;
			default:
				throw std::invalid_argument("Unknown start condition selected. Don't know how to generate proper observable file name.");
		}
	}
}

std::string meta::create_prng_name(const Inputparameters& parameters, int number) noexcept {
	using namespace std;
	std::stringstream strnumber;
	strnumber.fill('0');
	strnumber.width(parameters.get_config_number_digits());
	strnumber << right << number;
	stringstream outfilename;
	outfilename << parameters.get_prng_prefix() << strnumber.str() << parameters.get_prng_postfix();
	string outputfile = outfilename.str();
	return outputfile;
}

std::string meta::create_prng_name(const Inputparameters& parameters) noexcept {
	using namespace std;
	stringstream outfilename;
	outfilename << parameters.get_prng_prefix() << "save" << parameters.get_prng_postfix();
	string outputfile = outfilename.str();
	return outputfile;
}

std::string meta::create_profiling_data_filename(const Inputparameters& parameters, std::string executableName) noexcept {
	using namespace std;
	stringstream outfilename;
	outfilename << parameters.get_profiling_data_prefix() << executableName << parameters.get_profiling_data_postfix();
	string outputfile = outfilename.str();
	return outputfile;
}

std::pair<int, std::vector<const char*>> meta::addOptionsToArgv(int argc, const char** argv, std::vector<const char*> new_options)
{
	std::vector<const char*> argv_new(argv, argv + argc);
	argc += new_options.size();
	argv_new.insert(argv_new.end(), new_options.begin(), new_options.end());

	return std::make_pair(argc, argv_new);
}
