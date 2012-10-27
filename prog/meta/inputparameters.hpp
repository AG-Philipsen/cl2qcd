/** @file
 * Input file handling
 */

#ifndef _META_INPUTPARAMETERS_HPP_
#define _META_INPUTPARAMETERS_HPP_

#include <vector>
#include <string>

#include "../logger.hpp"

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
	class Inputparameters {

		public:

			enum action { wilson = 1, clover, twistedmass, tlsym, iwasaki, dbw2 };
			enum integrator { leapfrog = 1, twomn };
			enum startcondition { cold_start = 1, hot_start, start_from_source };
			enum solver { cg = 1, bicgstab, bicgstab_save };
	  enum sourcetype {point = 1, volume, timeslice};
	  enum sourcecontent {one = 1, z4, gaussian};

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

			/*
			 * Read accessor functions
			 */

			size_t get_precision() const noexcept;

			const std::vector<int> get_selected_devices() const noexcept;
			int get_device_count() const noexcept;
			bool get_use_gpu() const noexcept;
			bool get_use_cpu() const noexcept;

			bool get_use_aniso() const noexcept;
			bool get_use_chem_pot_re() const noexcept;
			bool get_use_chem_pot_im() const noexcept;
			bool get_use_smearing() const noexcept;
			bool get_use_mp() const noexcept;
			int get_nspace() const noexcept;
			int get_ntime() const noexcept;

			startcondition get_startcondition() const noexcept;
			bool get_saveconfigs() const noexcept;
			int get_writefrequency() const noexcept;
			int get_savefrequency() const noexcept;
			std::string get_sourcefile() const noexcept;
			bool get_print_to_screen() const noexcept;
			//This is obvious!!!
			uint32_t get_host_seed() const noexcept;

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
			double get_mu() const noexcept;
			double get_csw() const noexcept;
			int get_iter0() const noexcept;
			int get_iter1() const noexcept;
			double get_kappa_mp() const noexcept;
			double get_mu_mp() const noexcept;
			double get_csw_mp() const noexcept;
			int get_iter0_mp() const noexcept;
			int get_iter1_mp() const noexcept;
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
			bool get_use_pointsource() const noexcept;
			bool get_use_gauge_only() const noexcept;
			int get_num_sources() const noexcept;
			int get_pointsource_x() const noexcept;
			int get_pointsource_y() const noexcept;
			int get_pointsource_z() const noexcept;
			int get_pointsource_t() const noexcept;

			double get_solver_prec() const noexcept;
			double get_force_prec() const noexcept;
			int get_iter_refresh() const noexcept;
			int get_iter_refresh_mp() const noexcept;

			//HMC specific parameters
			double get_tau() const noexcept;
			bool get_reversibility_check() const noexcept;
			int get_integrationsteps(size_t timescale) const noexcept;
			int get_hmcsteps() const noexcept;
			int get_num_timescales() const noexcept;
			integrator get_integrator(size_t timescale) const noexcept;
			//this is the optimal value...
			double get_lambda(size_t timescale) const noexcept;

			//direction for the correlator
			int get_corr_dir() const noexcept;

			bool get_use_same_rnd_numbers() const noexcept;
			bool get_profile_solver() const noexcept;
			double get_test_ref_value() const noexcept;

			bool is_ocl_compiler_opt_disabled() const noexcept;

			bool get_use_merge_kernels_fermion() const noexcept;
			bool get_use_merge_kernels_spinor() const noexcept;
			bool get_use_rec12() const noexcept;

	                //parameters to read in gauge configurations
	  bool get_read_multiple_configs() const noexcept;
	  int get_config_read_start() const noexcept;
	  int get_config_read_end() const noexcept;
	  int get_config_read_incr() const noexcept;
	  int get_config_number_digits() const noexcept;
	  std::string get_config_prefix() const noexcept;
	  std::string get_config_postfix() const noexcept;
	  std::string get_ferm_obs_prefix() const noexcept;
	  std::string get_ferm_obs_postfix() const noexcept;
	  std::string get_gauge_obs_prefix() const noexcept;
	  std::string get_gauge_obs_postfix() const noexcept;
	  bool get_ferm_obs_to_single_file() const noexcept;
	  bool get_gauge_obs_to_single_file() const noexcept;
	  std::string get_hmc_obs_prefix() const noexcept;
	  std::string get_hmc_obs_postfix() const noexcept;
	  bool get_hmc_obs_to_single_file() const noexcept;

			std::string get_log_level() const noexcept;

		private:
			size_t precision;

			std::vector<int> selected_devices;
			int device_count;
			bool use_gpu;
			bool use_cpu;

			bool use_aniso;
			bool use_chem_pot_re;
			bool use_chem_pot_im;
			bool use_smearing;
			bool use_mp;
			int nspace;
			int ntime;

			startcondition _startcondition;
			bool saveconfigs;
			int writefrequency;
			int savefrequency;
			std::string sourcefile;
			bool print_to_screen;
			//This is obvious!!!
			uint32_t host_seed;

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
			double mu;
			double csw;
			int iter0;
			int iter1;
			double kappa_mp;
			double mu_mp;
			double csw_mp;
			int iter0_mp;
			int iter1_mp;
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
			bool use_pointsource;
			bool use_gauge_only;
			int num_sources;
			int pointsource_x;
			int pointsource_y;
			int pointsource_z;
			int pointsource_t;

			double solver_prec;
			double force_prec;
			int iter_refresh;
			int iter_refresh_mp;

			//HMC specific parameters
			double tau;
			bool reversibility_check;
			int integrationsteps0;
			int integrationsteps1;
			int integrationsteps2;
			int hmcsteps;
			int num_timescales;
			integrator integrator0;
			integrator integrator1;
			integrator integrator2;
			double lambda0;
			double lambda1;
			double lambda2;

			//direction for the correlator
			int corr_dir;

			bool use_same_rnd_numbers;
			bool profile_solver;
			double test_ref_value;

			bool ocl_compiler_opt_disabled;

			bool use_merge_kernels_fermion;
			bool use_merge_kernels_spinor;
			bool use_rec12;

	                //parameters to read in gauge configurations
	  bool read_multiple_configs;
	  int config_read_start;
	  int config_read_end;
	  int config_read_incr;
	  int config_number_digits;
	  std::string config_prefix;
	  std::string config_postfix;

	  //parameters to write out observables
	  bool gauge_obs_to_single_file;
	  std::string gauge_obs_prefix;
	  std::string gauge_obs_postfix;
	  bool ferm_obs_to_single_file;
	  std::string ferm_obs_prefix;
	  std::string ferm_obs_postfix;
	  bool hmc_obs_to_single_file;
	  std::string hmc_obs_prefix;
	  std::string hmc_obs_postfix;

			std::string log_level;
	};
}

#endif /* _META_INPUTPARAMETERS_H_ */
