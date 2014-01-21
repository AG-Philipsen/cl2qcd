/** @file
 * rhmc algorithm main()
 *
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

#include "rhmc.h"

#include "../meta/util.hpp"
#include <cmath>

//inserted from former general_header.h
#include <fstream>
using namespace std;

//these can be used to measure times on the host spent on gaugeobservables
usetimer poly_timer;
usetimer plaq_timer;

//timers for overall times
usetimer total_timer;
usetimer init_timer;
usetimer perform_timer;

//output-method for timer above and the copy-timer from (each) device
///@todo For hybrid case one will have to add more timer here
void general_time_output(usetimer * total, usetimer * init_timer, usetimer * perform_timer, usetimer * plaq_timer, usetimer * poly_timer);
// end of insertion

static void check_rhmc_parameters(const meta::Inputparameters& params);
static void print_rhmcobservables(const hmc_observables& obs, int iter, const std::string& filename, const meta::Inputparameters& params);
static void print_rhmcobservables(const hmc_observables& obs, int iter);

int main(int argc, const char* argv[])
{
	using physics::algorithms::perform_rhmc_step;
	using physics::algorithms::Rational_Approximation;

	try {
		//I build a vector from argv so that I can manually add parameters that I would
		//always add (e.g. fermact=rooted_stagg) without modifying default values in Inputparameters
		std::vector<const char*> w(argv, argv + argc);
		w.push_back("--fermact=rooted_stagg");
		const char ** _params = &w[0];
		meta::Inputparameters parameters(w.size(), _params);
		switchLogLevel(parameters.get_log_level());
		
		//Check if inputparameters has reasonable parameters
		check_rhmc_parameters(parameters);
		
		meta::print_info_rhmc(argv[0], parameters);

		ofstream ofile;
		ofile.open("rhmc.log");
		if(ofile.is_open()) {
			meta::print_info_rhmc(argv[0], &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for rhmc.";
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Initialization
		/////////////////////////////////////////////////////////////////////////////////
		init_timer.reset();
		hmc_observables obs;

		hardware::System system(parameters);
		physics::PRNG prng(system);

		logger.trace() << "Init gaugefield" ;
		physics::lattices::Gaugefield gaugefield(system, prng);

		logger.info() << "Gaugeobservables:";
		print_gaugeobservables(gaugefield, 0);
		init_timer.add();

		/////////////////////////////////////////////////////////////////////////////////
		// RHMC
		/////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		//start from the iterationnumber from sourcefile
		//NOTE: this is 0 in case of cold or hot start
		int iter = gaugefield.get_parameters_source().trajectorynr_source;
		const int rhmc_iter = iter + parameters.get_rhmcsteps();
		if(iter >= rhmc_iter)
		  logger.warn() << "The total number of RHMC iterations is NOT smaller than those already done.\nNO further iteration will be performed!";
		hmc_float acc_rate = 0.;
		const int writefreq = parameters.get_writefrequency();
		const int savefreq = parameters.get_savefrequency();
		
		logger.info() << "";
		logger.info() << "Generation of Rational Approximations...";
		Rational_Approximation *approx_hb, *approx_md, *approx_met;
		if(parameters.get_read_rational_approximations_from_file()){
			approx_hb  = new Rational_Approximation(parameters.get_approx_heatbath_file());
			approx_md  = new Rational_Approximation(parameters.get_approx_md_file());
			approx_met = new Rational_Approximation(parameters.get_approx_metropolis_file());
		}else{
			//This is the approx. to be used to generate the initial (pseudo)fermionic field
			approx_hb = new Rational_Approximation(parameters.get_metro_approx_ord(),
					parameters.get_num_tastes(), 8, parameters.get_approx_lower(),
					parameters.get_approx_upper(), false);
			//This is the approx. to be used to generate the initial (pseudo)fermionic field
			approx_md = new Rational_Approximation(parameters.get_md_approx_ord(),
					parameters.get_num_tastes(), 4, parameters.get_approx_lower(),
					parameters.get_approx_upper(), true);
			//This is the approx. to be used to generate the initial (pseudo)fermionic field
			approx_met = new Rational_Approximation(parameters.get_metro_approx_ord(),
					  parameters.get_num_tastes(), 4, parameters.get_approx_lower(),
					  parameters.get_approx_upper(), true);
		}
		
		
		logger.info() << "";
		logger.info() << "Perform RHMC on device(s)... ";

		//main hmc-loop
		for(; iter < rhmc_iter; iter++) {
			//generate new random-number for Metropolis step
			const hmc_float rnd_number = prng.get_double();
			
			obs = perform_rhmc_step(*approx_hb, *approx_md, *approx_met,
						 &gaugefield, iter, rnd_number, prng, system);

			acc_rate += obs.accept;
			if(((iter + 1)%writefreq) == 0) {
				std::string gaugeout_name = meta::get_rhmc_obs_file_name(parameters, "");
				print_rhmcobservables(obs, iter, gaugeout_name, parameters);
			} else if(parameters.get_print_to_screen()) {
				print_rhmcobservables(obs, iter);
			}

			if( savefreq != 0 && ( (iter + 1) % savefreq ) == 0 ) {
				// save gaugefield
				gaugefield.save(iter + 1);

				// save prng
				std::ostringstream prng_file_name;
				prng_file_name << "prng.";
				prng_file_name.fill('0');
				prng_file_name.width(parameters.get_config_number_digits());
				prng_file_name << iter + 1;
				prng.store(prng_file_name.str());
			}
		}

		//always save the config on the last iteration
		std::string outputfile = "conf.save";
		logger.info() << "saving current gaugefield to file \"" << outputfile << "\"";
		gaugefield.save(outputfile, iter);
		outputfile = "prng.save";
		logger.info() << "saving current prng state to \"" << outputfile << "\"";
		prng.store(outputfile);

		logger.info() << "RHMC done";
		logger.info() << "Acceptance rate: " << fixed <<  setprecision(1) << percent(acc_rate, parameters.get_rhmcsteps()) << "%";
		perform_timer.add();
		
		/////////////////////////////////////////////////////////////////////////////////
		// Final output
		/////////////////////////////////////////////////////////////////////////////////
		total_timer.add();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;

}


static void check_rhmc_parameters(const meta::Inputparameters& p)
{
	if(p.get_fermact() != meta::Inputparameters::rooted_stagg)
	  throw Invalid_Parameters("Fermion action not suitable for RHMC!", meta::Inputparameters::rooted_stagg, p.get_fermact());
	if(!p.get_use_eo())
	  throw Invalid_Parameters("RHMC available only WITH eo-prec!", "use_eo=1", p.get_use_eo());
	if(p.get_use_mp())
	  throw Invalid_Parameters("RHMC available only WITHOUT mass preconditionig!", "use_mp=0", p.get_use_mp());
	if(p.get_use_chem_pot_re())
	  throw Invalid_Parameters("RHMC available only WITHOUT real chemical potential!", "use_chem_pot_re=0", p.get_use_chem_pot_re());
	if(p.get_use_chem_pot_im())
	  throw Invalid_Parameters("RHMC available only WITHOUT imaginary chemical potential!", "use_chem_pot_im=0", p.get_use_chem_pot_im());
	if(p.get_num_tastes()%4 == 0)
	  throw Invalid_Parameters("RHMC not working with multiple of 4 tastes (there is no need of the rooting trick)!", "num_tastes%4 !=0", "num_tastes=" + to_string(p.get_num_tastes()));
	if(p.get_cg_iteration_block_size() == 0 || p.get_findminmax_iteration_block_size() == 0)
	  throw Invalid_Parameters("Iteration block sizes CANNOT be zero!", "cg_iteration_block_size!=0 && findminmax_iteration_block_size!=0", p.get_cg_iteration_block_size()==0 ? "cg_iteration_block_size=0" : "findminmax_iteration_block_size=0");
	if(p.get_approx_upper() != 1)
	  throw Invalid_Parameters("RHMC not available if the Rational expansion is not calculated in [..,1]!", "approx_upper=1", p.get_approx_upper());
	if(p.get_approx_lower() >= 1)
	   throw Invalid_Parameters("The lower bound of the Rational expansion >= than the upper one does not make sense!", "approx_lower < approx_upper", to_string(p.get_approx_lower())+" >= 1");
	if(p.get_writefrequency() == 0)
	  throw Invalid_Parameters("Write frequency CANNOT be zero!", "writefrequency!=0", p.get_writefrequency());
	if(p.get_savefrequency() == 0){
	  logger.warn() << "savefrequency==0 and hence NEITHER gaugefield NOR prng will be saved!!";
	  logger.warn() << "The simulation will start in 5 seconds..";
	  sleep(2); logger.warn() << "..3.."; sleep(1); logger.warn() << "..2.."; sleep(1);
	  logger.warn() << "..1, go!"; sleep(1);
	}
}





static void print_rhmcobservables(const hmc_observables& obs, const int iter, const std::string& filename, const meta::Inputparameters& params)
{
	const hmc_float exp_deltaH = std::exp(obs.deltaH);
	std::fstream hmcout(filename.c_str(), std::ios::out | std::ios::app);
	if(!hmcout.is_open()) throw File_Exception(filename);
	hmcout << iter << "\t";
	hmcout.width(8);
	hmcout.precision(15);
	//print plaquette (plaq, tplaq, splaq)
	hmcout << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq;
	//print polyakov loop (re, im, abs)
	hmcout << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" << sqrt(obs.poly.re * obs.poly.re + obs.poly.im * obs.poly.im);
	//print deltaH, exp(deltaH), acceptance-propability, accept (yes or no)
	hmcout <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept;
	//print number of iterations used in inversions with full and force precision
	/**
	 * @todo: The counters should be implemented once the solver class is used!"
	 * until then, only write "0"!
	 */
	int iter0 = 0;
	int iter1 = 0;
	hmcout << "\t" << iter0 << "\t" << iter1;
	if(params.get_use_mp() ) {
		hmcout << "\t" << iter0 << "\t" << iter1;
	}
	if(meta::get_use_rectangles(params) ) {
		//print rectangle value
		hmcout << "\t" << obs.rectangles;
	}
	hmcout << std::endl;
	hmcout.close();

	//print to screen
	print_rhmcobservables(obs, iter);
}

static void print_rhmcobservables(const hmc_observables& obs, const int iter)
{
	using namespace std;
	//short version of output, all obs are collected in the output file anyways...
	logger.info() << "\tRHMC [OBS]:\t" << iter << setw(8) << setfill(' ') << "\t" << setprecision(15) << obs.plaq << "\t" << obs.poly.re << "\t" << obs.poly.im;
}

void general_time_output(usetimer * total, usetimer * init_timer, usetimer * perform_timer, usetimer * plaq_timer, usetimer * poly_timer)
{

	uint64_t totaltime = (*total).getTime();

	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time

	uint64_t init_time = (*init_timer).getTime();
	uint64_t perform_time = (*perform_timer).getTime();
	uint64_t plaq_time = (*plaq_timer).getTime();
	uint64_t poly_time = (*poly_timer).getTime();

	int plaq_steps = (*plaq_timer).getNumMeas();
	int poly_steps = (*poly_timer).getNumMeas();

	uint64_t poly_avgtime = divide(poly_time, poly_steps);
	uint64_t plaq_avgtime = divide(plaq_time, plaq_steps);

	logger.info() << "## *******************************************************************";
	logger.info() << "## General Times [mus]:";
	logger.info() << "## *******************************************************************";
	logger.info() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
	logger.info() << "## Total:\t" << setfill(' ') << setw(12) << totaltime;
	logger.info() << "## Init.:\t" << setfill(' ') << setw(12) << init_time << '\t' << fixed << setw(5) << setprecision(1) << percent(init_time, totaltime) ;
	logger.info() << "## Perf.:\t" << setfill(' ') << setw(12) << perform_time << '\t' << fixed << setw(5) << setprecision(1) << percent(perform_time, totaltime) ;
	logger.info() << "## *******************************************************************";
	logger.info() << "## Host-Obs:\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg" << '\t' << setw(5) << "perc";
	logger.info() << "## Plaq.:\t" << setfill(' ') << setw(12) << plaq_time << '\t' << setw(12) << plaq_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(plaq_time, totaltime);
	logger.info() << "## Poly.:\t" << setfill(' ') << setw(12) << poly_time << '\t' << setw(12) << poly_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(poly_time, totaltime);
	logger.info() << "## *******************************************************************";

	logger.info() << "## writing general times to file: \"general_time_output\"";
	ofstream ofile;
	ofile.open("general_time_output");
	if(ofile.is_open()) {
		ofile  << "## *******************************************************************" << endl;
		ofile  << "## General Times [mus]:" << endl;
		ofile << "## Total\tInit\tPerformance\tHost-Plaq\tHost-Pol" << endl;
		ofile  << totaltime << "\t" << init_time << '\t' << perform_time << '\t' << plaq_time << '\t' << poly_time << endl;
		ofile.close();
	} else {
		logger.warn() << "Could not open output file for general time output.";
	}


	return;
}
