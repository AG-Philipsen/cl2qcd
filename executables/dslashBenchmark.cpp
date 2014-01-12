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
  spinorfield1 = new  hardware::buffers::Spinor(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device);
  spinorfield2 = new  hardware::buffers::Spinor(hardware::code::get_eoprec_spinorfieldsize(device->get_mem_lattice_size()), device);
}

void dslashBenchmark::performBenchmarkForSpecificKernels()
{
  auto gf_buffer = gaugefield->get_buffers().at(0);
  auto fermion_code = device->get_fermion_code();
  fermion_code->dslash_eo_device(spinorfield1, spinorfield2, gf_buffer, EVEN);
  fermion_code->dslash_eo_device(spinorfield1, spinorfield2, gf_buffer, ODD);
  device->synchronize();
}
