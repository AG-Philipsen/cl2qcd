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
#include "../physics/lattices/spinorfield_eo.hpp"
#include "../physics/fermionmatrix/fermionmatrix.hpp"
#include "../klepsydra/klepsydra.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/fermions.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
	using physics::fermionmatrix::dslash;

	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		//NOTE: This is inserted here because of the refactoring of the other executables
		logger.info() << "## Starting executable: " << argv[0];
		meta::print_info_global(parameters);
		meta::print_info_configs_io(parameters);
		meta::print_info_prng_io(parameters);
		meta::print_info_inverter(parameters);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		hardware::System system(parameters);

		physics::PRNG prng(system);
		physics::lattices::Gaugefield gf(system, prng);
		physics::lattices::Spinorfield_eo sf1(system);
		physics::lattices::Spinorfield_eo sf2(system);

		// update gaugefield buffers once to have update links fully initialized
		gf.update_halo();

		logger.info() << "Gaugeobservables:";
		print_gaugeobservables(gf, 0);

		dslash(&sf2, gf, sf1, EVEN);
		dslash(&sf1, gf, sf2, ODD);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// dslash-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for(auto dev: system.get_devices()) {
			dev->synchronize();
		}

		int hmc_iter = parameters.get_hmcsteps();

		logger.info() << "Perform dslash (EVEN + ODD) " << hmc_iter << " times.";
		klepsydra::Monotonic timer;
		for(int iter = 0; iter < hmc_iter; ++iter) {
			dslash(&sf2, gf, sf1, EVEN);
			dslash(&sf1, gf, sf2, ODD);
		}
		for(auto dev: system.get_devices()) {
			dev->synchronize();
		}
		auto elapsed_mus = timer.getTime();
		logger.trace() << "dslash benchmarking done" ;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		auto fermion_code = system.get_devices()[0]->get_fermion_code();
		size_t flop_count = fermion_code->get_flop_size("dslash_eo");
		size_t byte_count = fermion_code->get_read_write_size("dslash_eo");
		double gflops = static_cast<double>(flop_count) * 2 * hmc_iter / elapsed_mus / 1e3;
		double gbytes = static_cast<double>(byte_count) * 2 * hmc_iter / elapsed_mus / 1e3;
		logger.info() << "Dslash performance: " << gflops << " GFLOPS";
		logger.info() << "Dslash memory: " << gbytes << " GB/S";

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

