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
	void readfile(const char* ifn);

	/////////////////////////////////////////////////////
	// access to private members
	/**
	 * Reset all parameters to default values.
	 */
	void set_defaults();
	hmc_float get_kappa() const;
	void set_mubar_negative();
	void calc_mubar();
	hmc_float get_mubar() const;
	hmc_float get_beta() const;
	hmc_float get_tau() const;
	hmc_float get_theta_fermion_spatial() const;
	hmc_float get_theta_fermion_temporal() const;
	hmc_float get_theta_gaugefield() const;
	hmc_float get_mu() const;
	hmc_float get_csw() const;
	hmc_float get_chem_pot_re() const;
	hmc_float get_chem_pot_im() const;
	int get_nt() const;
	int get_ns() const;
	//	int get_xi() const;
	//	hmc_float get_xi_0() const;
	bool get_use_rec12() const;
	bool get_use_gpu() const;
	bool get_use_eo() const;
	int get_cgmax() const;
	/**
	 * The precision to be used for gaugefield storage in bits.
	 *
	 * @li 32 for single precision
	 * @li 64 for double precision
	 *
	 * @todo Currently this is also used in many places to specify the calculation
	 *       precision. However it should be possible to calculate in a different
	 *       precision than one uses for the input/output files.
	 */
	int get_prec() const;
	int get_startcondition() const;
	int get_thermalizationsteps() const;
	int get_heatbathsteps() const;
	int get_overrelaxsteps() const;
	int get_hmcsteps() const;
	int get_integrationsteps1() const;
	int get_integrationsteps2() const;
	hmc_float get_solver_prec() const;
	hmc_float get_force_prec() const;
	cl_uint get_iter_refresh() const;
	bool get_saveconfigs() const;
	int get_savefrequency() const;
	int get_writefrequency() const;
	int get_fermact() const;
	void display_sourcefile() const;
	void display_sourcefilenumber() const;
	int get_num_dev() const;
	/** Spatial volume of the lattice */
	int get_volspace() const;
	/** 4-Dimensional Volume of the lattice */
	int get_vol4d() const;
	int get_spinorsize() const;
	int get_halfspinorsize() const;
	int get_eoprec_spinorfieldsize() const;
	int get_gaugemomentasize() const;
	int get_su3algebrasize() const;
	int get_gaugefieldsize() const;
	int get_spinorfieldsize() const;
	/**
	 * Return the size of the buffer in bytes.
	 */
	int get_sf_buf_size() const;
	/**
	 * Return the size of the buffer in bytes.
	 */
	int get_eo_sf_buf_size() const;
	/**
	 * Return the size of the buffer in bytes.
	 */
	int get_gf_buf_size() const;
	/**
	 * Return the size of the buffer in bytes.
	 */
	int get_gm_buf_size() const;
	bool get_use_chem_pot_re() const;
	bool get_use_chem_pot_im() const;
	bool get_use_smearing() const;
	bool get_print_to_screen() const;
	int get_host_seed() const;
	hmc_float get_rho() const;
	int get_rho_iter() const;
	bool get_use_cg() const;
	bool get_use_bicgstab_save() const;
	bool get_use_autotuning() const;
	bool get_use_pointsource() const;
	bool get_reversibility_check() const;
	int get_num_sources() const;
	int get_source_pos_spatial() const;
	int get_source_pos_temporal() const;
	int get_corr_dir() const;
	int get_integrator() const;
	int get_num_timescales() const;
	hmc_float get_lambda1() const;
	hmc_float get_lambda2() const;

	bool get_use_same_rnd_numbers() const;
#ifdef _PROFILING_
	int get_mat_size() const;
	int get_float_size() const;
#endif
	/////////////////////////////////////////////////////
	// printing-functions for the different executables
	/**
	 * Print info for heatbath executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_heatbath(char* progname) const;

	/**
	 * Print info for heatbath executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_heatbath(char* progname, ostream* os) const;

	/**
	 * Print info for inverter executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_inverter(char* progname) const;

	/**
	 * Print info for inverter executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_inverter(char* progname, ostream* os) const;

	/**
	 * Print info for tk_kappa executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_tkkappa(char* progname) const;

	/**
	 * Print info for tk_kappa executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_tkkappa(char* progname, ostream* os) const;

	/**
	 * Print info for hmc executable using logger.
	 * @param progname Name of the executable.
	 */
	void print_info_hmc(char* progname) const;

	/**
	 * Print info for hmc executable using output stream.
	 * @param progname Name of the executable.
	 * @param os Pointer to wanted output stream.
	 */
	void print_info_hmc(char* progname, ostream* os) const;

	/**
	 * print info independent of executable using logger
	 */
	void print_info_global() const;

	/**
	 * print info independent of executable using output stream
	 */
	void print_info_global(ostream* os) const;

	/**
	 *
	 */
	void print_info_fermion() const;

	/**
	 *
	 */
	void print_info_fermion(ostream * os) const;

	/**
	 * check inputparameters against compile settings
	 * NOTE: In the end, this is propably not needed anymore, but for now it is a safety net
	 */
	void check_settings_global() const;

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
	//	int xi;
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
	hmc_float solver_prec;
	hmc_float force_prec;
	cl_uint iter_refresh;
	int rho_iter;
	bool reversibility_check;
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
	hmc_float lambda1;
	hmc_float lambda2;
	bool use_cg;
	bool use_bicgstab_save;
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
	void solver_assign(bool * out, std::string line);

	//this can bee used to generate rnd-fields with only one thread
	bool use_same_rnd_numbers;
};

#endif /* _INPUTPARAMETERSH_ */
