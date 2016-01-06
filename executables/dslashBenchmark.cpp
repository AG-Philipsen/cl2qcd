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
  spinorfield1 = new physics::lattices::Spinorfield_eo(*system, interfacesHandler->getInterface<physics::lattices::Spinorfield_eo>());
  spinorfield2 = new physics::lattices::Spinorfield_eo(*system, interfacesHandler->getInterface<physics::lattices::Spinorfield_eo>());
}

void dslashBenchmark::performBenchmarkForSpecificKernels()
{
  auto gaugefield_buffer = gaugefield->get_buffers().at(0);
  auto spinorfield1_buffer = spinorfield1->get_buffers().at(0);
  auto spinorfield2_buffer = spinorfield2->get_buffers().at(0);
  auto fermion_code = device->getFermionCode();
  fermion_code->dslash_eo_device(spinorfield1_buffer, spinorfield2_buffer, gaugefield_buffer, EVEN);
  fermion_code->dslash_eo_device(spinorfield1_buffer, spinorfield2_buffer, gaugefield_buffer, ODD);
  device->synchronize();
}

void dslashBenchmark::enqueueSpecificKernelForBenchmarkingMultipleDevices()
{
  physics::fermionmatrix::dslash(spinorfield2, *gaugefield, *spinorfield1, EVEN, parameters.get_kappa());
  physics::fermionmatrix::dslash(spinorfield1, *gaugefield, *spinorfield2, ODD, parameters.get_kappa());
}

void dslashBenchmark::printProfilingDataToScreen()
{
  auto fermion_code = system->get_devices()[0]->getFermionCode();
  size_t flop_count = fermion_code->get_flop_size("dslash_eo");
  size_t byte_count = fermion_code->get_read_write_size("dslash_eo");
  double gflops = static_cast<double>(flop_count) * 2 * benchmarkSteps / executionTime / 1e3;
  double gbytes = static_cast<double>(byte_count) * 2 * benchmarkSteps / executionTime / 1e3;
  logger.info() << "Dslash performance: " << gflops << " GFLOPS";
  logger.info() << "Dslash memory: " << gbytes << " GB/S";
}

