/*
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2018,2021 Alessandro Sciarra
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

#include "su3heatbathBenchmark.hpp"

#include "../hardware/code/heatbath.hpp"

su3heatbathBenchmark::su3heatbathBenchmark(int argc, const char* argv[])
    : benchmarkExecutable(argc, argv, "SU3_Heatbath"), overrelaxationSteps(1)
{
    if (system->get_devices().size() != 1) {
        logger.fatal() << "There must be exactly one device chosen for the heatbath benchmark to be performed.";
    }
    if (!parameters.get_enable_profiling()) {
        throw Print_Error_Message("Profiling is not enabled. Aborting...\n", __FILE__, __LINE__);
    }
}

void su3heatbathBenchmark::enqueueSpecificKernelForBenchmark()
{
    physics::algorithms::su3heatbath(*gaugefield, *prng, overrelaxationSteps);
}

std::vector<double> su3heatbathBenchmark::getExecutionTimesOnDevices()
{
    std::vector<double> times(devices.size(), 0.0);
    for (unsigned int i = 0; i < devices.size(); ++i) {
        auto code = devices[0]->getHeatbathCode();
        std::vector<cl_kernel> kernels{code->heatbath_even, code->heatbath_odd, code->overrelax_even,
                                       code->overrelax_odd};
        for (auto kernel : kernels)
            times[i] += devices[i]->getProfilingData(kernel).get_total_time();
    }
    return times;
}

size_t su3heatbathBenchmark::getFlopsPerKernelCall()
{
    auto code = devices[0]->getHeatbathCode();
    return code->get_flop_size("heatbath_even") + code->get_flop_size("heatbath_odd") +
           overrelaxationSteps * (code->get_flop_size("overrelax_even") + code->get_flop_size("overrelax_odd"));
    ;
}
size_t su3heatbathBenchmark::getMemoryPerKernelCall()
{
    auto code = devices[0]->getHeatbathCode();
    return code->get_read_write_size("heatbath_even") + code->get_read_write_size("heatbath_odd") +
           overrelaxationSteps *
               (code->get_read_write_size("overrelax_even") + code->get_read_write_size("overrelax_odd"));
}
