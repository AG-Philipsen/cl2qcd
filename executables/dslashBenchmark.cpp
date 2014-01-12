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

/*
 * @file
 * Declaration of the hmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Hybrid Monte Carlo (HMC) algorithm.
 */

#include "dslashBenchmark.h"

dslashBenchmark::dslashBenchmark(int argc, const char* argv[]) :
  benchmarkExecutable(argc, argv)
{
  if(system->get_devices().size() != 1) {
    logger.fatal() << "There must be exactly one device chosen for the dslash benchmark to be performed.";
  }
  if(! parameters.get_enable_profiling() )
    {
      throw Print_Error_Message( "Profiling is not enabled. Aborting...\n", __FILE__, __LINE__);
    }
  spinorfield1 = new physics::lattices::Spinorfield_eo(*system);
  spinorfield2 = new physics::lattices::Spinorfield_eo(*system);
}

void dslashBenchmark::performBenchmarkForSpecificKernels()
{
  auto gaugefield_buffer = gaugefield->get_buffers().at(0);
  auto spinorfield1_buffer = spinorfield1->get_buffers().at(0);
  auto spinorfield2_buffer = spinorfield2->get_buffers().at(0);
  auto fermion_code = device->get_fermion_code();
  fermion_code->dslash_eo_device(spinorfield1_buffer, spinorfield2_buffer, gaugefield_buffer, EVEN);
  fermion_code->dslash_eo_device(spinorfield1_buffer, spinorfield2_buffer, gaugefield_buffer, ODD);
  device->synchronize();
}

void dslashBenchmark::benchmarkMultipleDevices()
{
  
		// update gaugefield buffers once to have update links fully initialized
		gaugefield->update_halo();

		physics::fermionmatrix::dslash(spinorfield2, *gaugefield, *spinorfield1, EVEN);
		physics::fermionmatrix::dslash(spinorfield1, *gaugefield, *spinorfield2, ODD);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// dslash-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for(auto dev: system->get_devices()) {
			dev->synchronize();
		}

		logger.info() << "Perform dslash (EVEN + ODD) " << benchmarkSteps << " times.";
		klepsydra::Monotonic timer;
		for(int iteration = 0; iteration < benchmarkSteps; ++iteration) {
		  physics::fermionmatrix::dslash(spinorfield2, *gaugefield, *spinorfield1, EVEN);
		  physics::fermionmatrix::dslash(spinorfield1, *gaugefield, *spinorfield2, ODD);
		}
		for(auto dev: system->get_devices()) {
			dev->synchronize();
		}
		auto elapsed_mus = timer.getTime();
		logger.trace() << "dslash benchmarking done" ;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		auto fermion_code = system->get_devices()[0]->get_fermion_code();
		size_t flop_count = fermion_code->get_flop_size("dslash_eo");
		size_t byte_count = fermion_code->get_read_write_size("dslash_eo");
		double gflops = static_cast<double>(flop_count) * 2 * benchmarkSteps / elapsed_mus / 1e3;
		double gbytes = static_cast<double>(byte_count) * 2 * benchmarkSteps / elapsed_mus / 1e3;
		logger.info() << "Dslash performance: " << gflops << " GFLOPS";
		logger.info() << "Dslash memory: " << gbytes << " GB/S";
  
}
