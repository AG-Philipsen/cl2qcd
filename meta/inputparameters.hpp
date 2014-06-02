/** @file
 * Input file handling
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

#ifndef _META_INPUTPARAMETERS_HPP_
#define _META_INPUTPARAMETERS_HPP_

#include "../host_functionality/logger.hpp"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;

#include "parametersIo.hpp"
#include "parametersObs.hpp"
#include "parametersHmc.hpp"
#include "parametersRhmc.hpp"
#include "parametersTest.hpp"

/**
 * This namespace contains generic utility code required by the other packages.
 */
namespace meta {


/**
 * Parser and representation of an input file.
 *
 * This class is copyable and assignable, but should
 * be used as a const value after initialization.
 */
	class Inputparameters : 
		public ParametersIo,
		public ParametersObs, 
		public ParametersHmc, 
		public ParametersRhmc,
		public ParametersTest
	{

public:

	enum action { wilson = 1, clover, twistedmass, tlsym, iwasaki, dbw2, rooted_stagg };
	enum startcondition { cold_start = 1, hot_start, start_from_source };
	enum solver { cg = 1, bicgstab, bicgstab_save };
	enum sourcetypes {point = 1, volume, timeslice, zslice};
	enum sourcecontents {one = 1, z4, gaussian, z2};

	/**
	 * The parsing of the input parameters aborted for some reason.
	 * Could either be invalid input or the specification of --help.
	 */
	struct parse_aborted {};

	/**
	 * Construct from command line and config file.
	 *
	 * Config file will be retrieved from command line.
	 *
	 * @throws parse_aborted
	 */
	Inputparameters(int argc, const char** argv);

	po::options_description getParameters_hmc();

	/*
	 * Read accessor functions
	 */

	size_t get_precision() const noexcept;

	const std::vector<int> get_selected_devices() const noexcept;
	int get_device_count() const noexcept;
	bool get_use_gpu() const noexcept;
	bool get_use_cpu() const noexcept;
	bool get_enable_profiling() const noexcept;

	bool get_use_aniso() const noexcept;
	bool get_use_chem_pot_re() const noexcept;
	bool get_use_chem_pot_im() const noexcept;
	bool get_use_smearing() const noexcept;
	int get_nspace() const noexcept;
	int get_ntime() const noexcept;

	int get_benchmarksteps() const noexcept;

	startcondition get_startcondition() const noexcept;
	std::string get_sourcefile() const noexcept;
	bool get_ignore_checksum_errors() const noexcept;
	bool get_print_to_screen() const noexcept;
	//This is obvious!!!
	uint32_t get_host_seed() const noexcept;
	std::string get_initial_prng_state() const noexcept;

	//gaugefield parameters
	double get_beta() const noexcept;
	double get_rho() const noexcept;
	int get_rho_iter() const noexcept;
	action get_gaugeact() const noexcept;

	//heatbath parameters
	int get_thermalizationsteps() const noexcept;
	int get_heatbathsteps() const noexcept;
	int get_overrelaxsteps() const noexcept;
	int get_xi() const noexcept;


	//fermionic parameters
	action get_fermact() const noexcept;
	action get_fermact_mp() const noexcept;
	double get_kappa() const noexcept;
	double get_mass() const noexcept;
	double get_mu() const noexcept;
	double get_csw() const noexcept;
	double get_kappa_mp() const noexcept;
	double get_mu_mp() const noexcept;
	double get_csw_mp() const noexcept;
	int get_cgmax() const noexcept;
	int get_cgmax_mp() const noexcept;
	double get_theta_fermion_spatial() const noexcept;
	double get_theta_fermion_temporal() const noexcept;
	double get_chem_pot_re() const noexcept;
	double get_chem_pot_im() const noexcept;
	bool get_use_eo() const noexcept;
	//at the moment, only 2 solvers are implemented..
	solver get_solver() const noexcept;
	solver get_solver_mp() const noexcept;
	int get_num_sources() const noexcept;
	int get_source_x() const noexcept;
	int get_source_y() const noexcept;
	int get_source_z() const noexcept;
	int get_source_t() const noexcept;
	bool get_place_sources_on_host() const noexcept;

	double get_solver_prec() const noexcept;
	double get_force_prec() const noexcept;
	int get_iter_refresh() const noexcept;
	int get_iter_refresh_mp() const noexcept;

	bool get_use_same_rnd_numbers() const noexcept;
	bool get_profile_solver() const noexcept;


	bool is_ocl_compiler_opt_disabled() const noexcept;

	bool get_use_merge_kernels_fermion() const noexcept;
	bool get_use_merge_kernels_spinor() const noexcept;
	bool get_use_rec12() const noexcept;

	//parameters to read in gauge configurations
	bool get_read_multiple_configs() const noexcept;
	int get_config_read_start() const noexcept;
	int get_config_read_end() const noexcept;
	int get_config_read_incr() const noexcept;

	std::string get_log_level() const noexcept;

	sourcetypes get_sourcetype() const noexcept;
	sourcecontents get_sourcecontent() const noexcept;


	int get_cg_iteration_block_size() const noexcept;
	bool get_cg_use_async_copy() const noexcept;
	int get_cg_minimum_iteration_count() const noexcept;

	bool get_split_cpu() const noexcept;

private:
	size_t precision;

	std::vector<int> selected_devices;
	int device_count;
	bool use_gpu;
        bool use_cpu;
        bool enable_profiling;

	bool use_aniso;
	bool use_chem_pot_re;
	bool use_chem_pot_im;
	bool use_smearing;
	int nspace;
	int ntime;

	startcondition _startcondition;

	std::string sourcefile;
	bool ignore_checksum_errors;
	bool print_to_screen;
	//This is obvious!!!
	uint32_t host_seed;
	std::string initial_prng_state;

	//gaugefield parameters
	double beta;
	double rho;
	int rho_iter;
	action gaugeact;

	//heatbath parameters
	int thermalizationsteps;
	int heatbathsteps;
	int overrelaxsteps;
	int xi;

	//fermionic parameters
	action fermact;
	action fermact_mp;
	double kappa;
	double mass; //staggered quark mass
	double mu;
	double csw;
	double kappa_mp;
	double mu_mp;
	double csw_mp;
	int cgmax;
	int cgmax_mp;
	double theta_fermion_spatial;
	double theta_fermion_temporal;
	double chem_pot_re;
	double chem_pot_im;
	bool use_eo;
	//at the moment, only 2 solvers are implemented..
	solver _solver;
	solver _solver_mp;
	int num_sources;
	int source_x;
	int source_y;
	int source_z;
	int source_t;
	bool place_sources_on_host;

	double solver_prec;
	double force_prec;
	int iter_refresh;
	int iter_refresh_mp;

	bool use_same_rnd_numbers;
	bool profile_solver;


	bool ocl_compiler_opt_disabled;

	bool use_merge_kernels_fermion;
	bool use_merge_kernels_spinor;
	bool use_rec12;

	//parameters to read in gauge configurations
	bool read_multiple_configs;
	int config_read_start;
	int config_read_end;
	int config_read_incr;

	std::string log_level;

	sourcetypes sourcetype;
	sourcecontents sourcecontent;

	int cg_iteration_block_size;
	bool cg_use_async_copy;
	int cg_minimum_iteration_count;

	bool split_cpu;
};
}

#endif /* _META_INPUTPARAMETERS_H_ */
