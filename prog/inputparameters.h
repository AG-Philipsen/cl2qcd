/** @file
 * Input file handling
 */

#ifndef _INPUTPARAMETERSH_
#define _INPUTPARAMETERSH_

//CP: this should know all types that might be used in the programm
#include "types.h"
#include "types_fermions.h"
#include "types_hmc.h"

#include "globaldefs.h"
#include "logger.hpp"

#include "exceptions.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#define EXIT_INPUTPARAMETERS 1

#include "host_geometry.h"

using namespace std;



/**
 * Parser for the input file.
 */
class inputparameters {
public:
	/**
	 * Default constructor loads default values.
	 */
	inputparameters() {
		set_defaults();
	};
	/**
	 * Parse the given file, overwriting current
	 * parameter values.
	 *
	 * \param ifn Name of the input file to parse
	 */
	void readfile(char* ifn);

	/////////////////////////////////////////////////////
	// access to private members
	/**
	 * Reset all parameters to default values.
	 */
	void set_defaults();
	hmc_float get_kappa();
	void set_mubar_negative();
	void calc_mubar();
	hmc_float get_mubar();
	hmc_float get_beta();
	hmc_float get_tau();
	hmc_float get_theta_fermion_spatial();
	hmc_float get_theta_fermion_temporal();
	hmc_float get_theta_gaugefield();
	hmc_float get_mu();
	hmc_float get_csw();
	hmc_float get_chem_pot_re();
	hmc_float get_chem_pot_im();
	int get_nt();
	int get_ns();
	bool get_use_rec12();
	bool get_use_gpu();
	bool get_use_eo();
	int get_cgmax();
	/**
	 * The precision to be used for gaugefield storage in bits.
	 * @li 32 for single precision
	 * @li 64 for double precision
	 */
	int get_prec();
	int get_startcondition();
	int get_thermalizationsteps();
	int get_heatbathsteps();
	int get_overrelaxsteps();
	int get_hmcsteps();
	int get_integrationsteps1();
	int get_integrationsteps2();
	bool get_saveconfigs();
	int get_savefrequency();
	int get_writefrequency();
	int get_fermact();
	void display_sourcefile();
	void display_sourcefilenumber();
	int get_num_dev();
	int get_volspace();
	int get_vol4d();
	int get_spinorsize();
	int get_halfspinorsize();
	int get_eoprec_spinorfieldsize();
	int get_gaugemomentasize();
	int get_su3algebrasize();
	int get_gaugefieldsize();
	int get_spinorfieldsize();
	int get_sf_buf_size();
	int get_eo_sf_buf_size();
	int get_gf_buf_size();
	int get_gm_buf_size();
	bool get_use_chem_pot_re();
	bool get_use_chem_pot_im();
	bool get_use_smearing();
	bool get_print_to_screen();
	int get_host_seed();
	hmc_float get_rho();
	int get_rho_iter();
	bool get_use_cg();
	bool get_use_autotuning();
	bool get_use_pointsource();
	int get_num_sources();
	int get_source_pos_spatial();
	int get_source_pos_temporal();
	int get_corr_dir();
	int get_integrator();
	int get_num_timescales();
#ifdef _PROFILING_
	int get_mat_size();
	int get_float_size();
#endif
	/////////////////////////////////////////////////////
	// printing-functions for the different executables
	/**
	 * Print info for heatbath executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_heatbath(char* progname);

	/**
	 * Print info for heatbath executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_heatbath(char* progname, ostream* os);

	/**
	 * Print info for inverter executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_inverter(char* progname);

	/**
	 * Print info for inverter executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_inverter(char* progname, ostream* os);

	/**
	 * Print info for tk_kappa executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_tkkappa(char* progname);

	/**
	 * Print info for tk_kappa executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_tkkappa(char* progname, ostream* os);

	/**
	 * Print info for hmc executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_hmc(char* progname);

	/**
	 * Print info for hmc executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_hmc(char* progname, ostream* os);

	/**
	 * print info independent of executable using logger
	 */
	void print_info_global();

	/**
	 * print info independent of executable using output stream
	 */
	void print_info_global(ostream* os);

	/**
	 *
	 */
	void print_info_fermion();

	/**
	 *
	 */
	void print_info_fermion(ostream * os);

	/**
	 * check inputparameters against compile settings
	 * NOTE: In the end, this is propably not needed anymore, but for now it is a safety net
	 */
	void check_settings_global();

	/**
	 * set global settings according to inputparameters
	 */
	void set_settings_global();

	//CP
	//this is out of laziness
	std::string sourcefile;
	std::string sourcefilenumber;
private:
	//general parameters
	int nspace;
	int ntime;
	int volspace;
	int vol4d;
	int spinorsize;
	int halfspinorsize;
	int spinorfieldsize;
	int eoprec_spinorfieldsize;
	int gaugemomentasize;
	int su3algebrasize;
	int gaugefieldsize;
	//sizes of fieldbuffers
	int sf_buf_size;
	int eo_sf_buf_size;
	int gf_buf_size;
	int gm_buf_size;
	
	bool use_rec12;
	bool use_gpu;
	bool use_eo;
	bool use_chem_pot_re;
	bool use_chem_pot_im;
	bool use_smearing;
	bool print_to_screen;
	bool use_autotuning;

#ifdef _PROFILING_
	//parameters that describe the size of datatypes in bytes
	int mat_size;
	int float_size;
#endif

	//more specific ones
	hmc_float kappa;
	hmc_float beta;
	hmc_float mu;
	hmc_float mubar;
	hmc_float csw;
	hmc_float theta_fermion_spatial;
	hmc_float theta_fermion_temporal;
	hmc_float theta_gaugefield;
	hmc_float chem_pot_re;
	hmc_float chem_pot_im;
	hmc_float tau;
	hmc_float rho;
	int rho_iter;
	long long host_seed;
	int num_dev;
	int cgmax;
	int fermact;
	int prec;
	int startcondition;
	int thermalizationsteps;
	int heatbathsteps;
	int overrelaxsteps;
	int hmcsteps;
	int savefrequency;
	bool saveconfigs;
	int writefrequency;
	int num_timescales;
	int integrationsteps1;
	int integrationsteps2;
	bool use_cg;
	bool use_pointsource;
	int num_sources;
	int pointsource_x;
	int pointsource_y;
	int pointsource_z;
	int pointsource_t;
	int corr_dir;
	int integrator;
	void val_assign(hmc_float* out, std::string line);
	void val_assign(int * out, std::string line);
	void sourcefilenumber_assign(std::string * out);
	void startcond_assign(int * out, std::string line);
	void fermact_assign(int * out, std::string line);
	void integrator_assign(int * out, std::string line);
	void val_assign(std::string * out, std::string line);
	void bool_assign(bool * out, std::string line);
};

#endif /* _INPUTPARAMETERSH_ */
