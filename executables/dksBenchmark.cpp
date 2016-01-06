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

#include "dksBenchmark.h"

dksBenchmark::dksBenchmark(int argc, const char* argv[]) : benchmarkExecutable(argc, argv)
{
  staggeredfield1 = new physics::lattices::Staggeredfield_eo(*system, interfacesHandler->getInterface<physics::lattices::Staggeredfield_eo>());
  staggeredfield2 = new physics::lattices::Staggeredfield_eo(*system, interfacesHandler->getInterface<physics::lattices::Staggeredfield_eo>());
}

void dksBenchmark::performBenchmarkForSpecificKernels()
{
  auto gaugefield_buffer = gaugefield->get_buffers().at(0);
  auto staggeredfield1_buffer = staggeredfield1->get_buffers().at(0);
  auto staggeredfield2_buffer = staggeredfield2->get_buffers().at(0);
  auto fermion_staggered_code = device->getFermionStaggeredCode();
  fermion_staggered_code->D_KS_eo_device(staggeredfield1_buffer, staggeredfield2_buffer, gaugefield_buffer, EVEN);
  fermion_staggered_code->D_KS_eo_device(staggeredfield1_buffer, staggeredfield2_buffer, gaugefield_buffer, ODD);
  device->synchronize();
}

void dksBenchmark::enqueueSpecificKernelForBenchmarkingMultipleDevices()
{
  physics::fermionmatrix::DKS_eo(staggeredfield2, *gaugefield, *staggeredfield2, EVEN);
  physics::fermionmatrix::DKS_eo(staggeredfield2, *gaugefield, *staggeredfield2, ODD);
}

void dksBenchmark::printProfilingDataToScreen()
{
  auto fermion_staggered_code = system->get_devices()[0]->getFermionStaggeredCode();
  size_t flop_count = fermion_staggered_code->get_flop_size("D_KS_eo");
  size_t byte_count = fermion_staggered_code->get_read_write_size("D_KS_eo");
  double gflops = static_cast<double>(flop_count) * 2 * benchmarkSteps / executionTime / 1e3;
  double gbytes = static_cast<double>(byte_count) * 2 * benchmarkSteps / executionTime / 1e3;
  logger.info() << "D_KS_eo performance: " << gflops << " GFLOPS";
  logger.info() << "D_KS_eo memory: " << gbytes << " GB/S";
  logger.info() << "Measured TIME: " << executionTime / 1.e3 << "msec";
}

