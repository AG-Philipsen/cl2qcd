/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "../hmc.h"

#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/algorithms/hmc.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
	meta::Inputparameters parameters(argc, argv);
	switchLogLevel(parameters.get_log_level());

	meta::print_info_hmc(argv[0], parameters);

	//name of file to store gauge observables, print initial information
	/** @todo think about what is a senseful filename*/
	stringstream gaugeout_name;
	gaugeout_name << "hmc_output";

	fstream logfile;
	logfile.open("hmc.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
		meta::print_info_hmc(argv[0], &logfile, parameters);
		logfile.close();
	} else {
		logger.warn() << "Could not open hmc.log";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//cl_int err;

	init_timer.reset();
	sourcefileparameters parameters_source;
	//hmc_observables obs;

	hardware::System system(parameters, true);
	physics::PRNG prng(system);

	logger.trace() << "init gaugefield" ;
	physics::lattices::Gaugefield gaugefield(system, prng);

	logger.trace() << "initial gaugeobservables:";
	print_gaugeobservables(gaugefield, 0);

	init_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HMC-Benchmarks
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	perform_timer.reset();
	//CP: this is taken from gaugefield_hmc. However, I took out all parts out of the loop that have to do with changing the fields in order to always perform the same HMC step.
	const int hmc_iter = parameters.get_hmcsteps();
	int iter;
	//This is the random-number generator for the metropolis-step
	const hmc_float rnd_number = prng.get_double();
	hmc_observables obs;

	throw Print_Error_Message("Not implemented");

//	logger.debug() << "\tinit spinorfield and gaugemomentum" ;
//	gaugefield.init_gaugemomentum_spinorfield(&solver_timer, prng);
//
//	logger.debug() << "\tupdate gaugefield and gaugemomentum" ;
//	//copy u->u' p->p' for the integrator
//	hardware::buffers::copyData(gaugefield.get_task_hmc(0)->get_new_u(), gaugefield.get_task_hmc(0)->get_device()->get_gaugefield_code()->get_gaugefield());
//	hardware::buffers::copyData(gaugefield.get_task_hmc(0)->get_clmem_new_p(), gaugefield.get_task_hmc(0)->get_clmem_p());
//	logger.trace() << "Perform " << hmc_iter << "of benchmarking";
//	for(iter = 0; iter < hmc_iter; iter ++) {
//		//here, clmem_phi is inverted several times and stored in clmem_phi_inv
//		gaugefield.integrator(&solver_timer);
//		//metropolis step: afterwards, the updated config is again in gaugefield and p
//		logger.debug() << "\tperform Metropolis step: " ;
//		//this call calculates also the HMC-Observables
//		obs = gaugefield.get_task_hmc(0)->metropolis(rnd_number, gaugefield.get_parameters().get_beta(), gaugefield.get_task_hmc(0)->get_device()->get_gaugefield_code()->get_gaugefield());
//		//CP: just reject the outcome of the metropolis step
//		logger.trace() << "\tfinished HMC trajectory " << iter ;
//	}
//	logger.trace() << "HMC-benchmarking done";
//	perform_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

	//CP: this is just a fist version and will go into an own file later
	string profiling_out;
	profiling_out = string(argv[0]) + string("_profiling_data");

	fstream prof_file;
	prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
		meta::print_info_heatbath(argv[0], &prof_file, parameters);
		prof_file.close();
	} else {
		logger.warn() << "Could not open " << profiling_out;
	}
	print_profiling(system, profiling_out);

	return 0;
}
