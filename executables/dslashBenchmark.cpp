/*
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2015,2018,2021 Alessandro Sciarra
 * Copyright (c) 2015 Christopher Czaban
 * Copyright (c) 2015 Francesca Cuteri
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @file
 * Declaration of the hmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Hybrid Monte Carlo (HMC) algorithm.
 */

#include "dslashBenchmark.hpp"

dslashBenchmark::dslashBenchmark(int argc, const char* argv[]) : benchmarkExecutable(argc, argv, "Dslash")
{
    spinorfield1       = new physics::lattices::Spinorfield_eo(*system,
                                                         interfacesHandler
                                                             ->getInterface<physics::lattices::Spinorfield_eo>());
    spinorfield2       = new physics::lattices::Spinorfield_eo(*system,
                                                         interfacesHandler
                                                             ->getInterface<physics::lattices::Spinorfield_eo>());
    parameters.fermact = common::action::wilson;
}

void dslashBenchmark::enqueueSpecificKernelForBenchmark()
{
    physics::fermionmatrix::dslash(spinorfield2, *gaugefield, *spinorfield1, EVEN, parameters.get_kappa());
    physics::fermionmatrix::dslash(spinorfield1, *gaugefield, *spinorfield2, ODD, parameters.get_kappa());
}

std::vector<double> dslashBenchmark::getExecutionTimesOnDevices()
{
    std::vector<double> times(devices.size(), std::numeric_limits<double>::quiet_NaN());
    for (unsigned int i = 0; i < devices.size(); ++i) {
#ifdef ASYNC_HALO_UPDATES
        auto code = devices[i]->getFermionCode();
        std::vector<cl_kernel> kernels{code->_dslash_eo_inner, code->_dslash_eo_boundary};
        times[i] = 0.0;
        for (auto kernel : kernels)
            times[i] += devices[i]->getProfilingData(kernel).get_total_time();
#else
        auto kernel = devices[i]->getFermionCode()->dslash_eo;
        times[i]    = devices[i]->getProfilingData(kernel).get_total_time();
#endif
    }
    return times;
}

size_t dslashBenchmark::getFlopsPerKernelCall()
{
    return devices[0]->getFermionCode()->get_flop_size("dslash_eo");
}

size_t dslashBenchmark::getMemoryPerKernelCall()
{
    return devices[0]->getFermionCode()->get_read_write_size("dslash_eo");
}
