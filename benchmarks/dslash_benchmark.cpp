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

#include <fstream>
#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/fermions.hpp"
#include "../hardware/code/spinors.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
	using namespace hardware::buffers;

	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		meta::print_info_inverter(argv[0], parameters);

		if(parameters.get_profile_solver() == false) {
			logger.warn() << "solver times will not be measured!";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		hardware::System system(parameters, true);

		// expect there to be exactly one device
		if(system.get_devices().size() != 1) {
			logger.fatal() << "There must be exactly one device chosen for the dslash benchmark to be performed.";
		}

		physics::PRNG prng(system);
		physics::lattices::Gaugefield gaugefield(system, prng);

		logger.info() << "Gaugeobservables:";
		print_gaugeobservables(gaugefield, 0);

		//init needed buffers again
		//these are: 2 eoprec spinorfield, 1 gaugefield
		auto gf_buffer = gaugefield.get_buffers().at(0);
		auto device = gf_buffer->get_device();
		const Spinor sf1(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device);
		const Spinor sf2(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// dslash-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		int hmc_iter = parameters.get_hmcsteps();
		int iter;
		auto fermion_code = device->get_fermion_code();

		logger.trace() << "Perform " << hmc_iter << "of dslash benchmarking (EVEN + ODD) for each step";
		for(iter = 0; iter < hmc_iter; iter ++) {
			fermion_code->dslash_eo_device(&sf1, &sf2, gf_buffer, EVEN);
			fermion_code->dslash_eo_device(&sf1, &sf2, gf_buffer, ODD);
		}
		device->synchronize();
		logger.trace() << "dslash benchmarking done" ;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//CP: this is just a fist version and will go into an own file later
		const std::string profiling_out = std::string(argv[0]) + std::string("_profiling_data");

		std::fstream prof_file;
		prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
		if(prof_file.is_open()) {
			meta::print_info_inverter(argv[0], &prof_file, parameters);
			prof_file.close();
		} else {
			logger.warn() << "Could not open " << profiling_out;
		}

		fermion_code->print_profiling(profiling_out, 0) ;

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
	}
}
