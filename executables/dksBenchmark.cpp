/*
 * Copyright (c) 2014,2015,2018,2021 Alessandro Sciarra
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

#include "dksBenchmark.hpp"

#include <numeric>

dksBenchmark::dksBenchmark(int argc, const char* argv[]) : benchmarkExecutable(argc, argv, "D_KS")
{
    staggeredfield1 = std::make_unique<
        physics::lattices::Staggeredfield_eo>(*system,
                                              interfacesHandler->getInterface<physics::lattices::Staggeredfield_eo>());
    staggeredfield2 = std::make_unique<
        physics::lattices::Staggeredfield_eo>(*system,
                                              interfacesHandler->getInterface<physics::lattices::Staggeredfield_eo>());
    parameters.fermact = common::action::rooted_stagg;
}

void dksBenchmark::enqueueSpecificKernelForBenchmark()
{
    physics::fermionmatrix::DKS_eo(staggeredfield2.get(), *gaugefield, *staggeredfield1, EVEN);
    physics::fermionmatrix::DKS_eo(staggeredfield2.get(), *gaugefield, *staggeredfield1, ODD);
}

std::vector<double> dksBenchmark::getExecutionTimesOnDevices()
{
    std::vector<double> times(devices.size(), std::numeric_limits<double>::quiet_NaN());
    for (unsigned int i = 0; i < devices.size(); ++i) {
        auto kernel = devices[i]->getFermionStaggeredCode()->D_KS_eo;
        times[i]    = devices[i]->getProfilingData(kernel).get_total_time();
    }
    return times;
}

size_t dksBenchmark::getFlopsPerKernelCall()
{
    return devices[0]->getFermionStaggeredCode()->get_flop_size("D_KS_eo");
}

size_t dksBenchmark::getMemoryPerKernelCall()
{
    return devices[0]->getFermionStaggeredCode()->get_read_write_size("D_KS_eo");
}
