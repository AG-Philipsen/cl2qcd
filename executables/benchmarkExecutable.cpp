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

#include "benchmarkExecutable.h"

benchmarkExecutable::benchmarkExecutable(int argc, const char* argv[]) : generalExecutable (argc, argv)
{
	initializationTimer.reset();
	gaugefield = new physics::lattices::Gaugefield(*system, &(interfacesHandler->getInterface<physics::lattices::Gaugefield>()), *prng);
	device = system->get_devices().at(0);
	benchmarkSteps = parameters.get_benchmarksteps();
	executionTime = 0;
	initializationTimer.add();
}

void benchmarkExecutable::benchmark()
{
    if(! parameters.get_enable_profiling() )
    {
      throw Print_Error_Message( "Profiling is not enabled. Aborting...\n", __FILE__, __LINE__);
    }
    if(system->get_devices().size() != 1) 
      {
	throw Print_Error_Message("There must be exactly one device chosen for this benchmark to be performed. Aborting...\n", __FILE__, __LINE__);
      }
    logger.info() << "Perform benchmarks..";
    performanceTimer.reset();
    for (int iteration = 0; iteration < benchmarkSteps; iteration++)
      {
	performBenchmarkForSpecificKernels();
      }
    performanceTimer.add();
    logger.info() << "Benchmarks done";
}

void benchmarkExecutable::benchmarkMultipleDevices()
{
    performanceTimer.reset();
    if( parameters.get_enable_profiling() )
    {
      throw Print_Error_Message( "Profiling is enabled. Aborting...\n", __FILE__, __LINE__);
    }
  // update gaugefield buffers once to have update links fully initialized
  gaugefield->update_halo();
  // ensure that the kernels are already built
  enqueueSpecificKernelForBenchmarkingMultipleDevices();
  
  synchronizeAllDevices();
  
  logger.info() << "Perform " << benchmarkSteps << " benchmarking steps.";
  klepsydra::Monotonic timer;
  for(int iteration = 0; iteration < benchmarkSteps; ++iteration) {
    enqueueSpecificKernelForBenchmarkingMultipleDevices();
  }
  synchronizeAllDevices();
  
  executionTime = timer.getTime();
  logger.info() << "Benchmarking done" ;
  
  printProfilingDataToScreen();
  performanceTimer.add();
}

void benchmarkExecutable::synchronizeAllDevices()
{
  for(auto dev: system->get_devices()) {
    dev->synchronize();
  }
}

