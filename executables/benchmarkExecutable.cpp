/*
 * Copyright (c) 2014-2016,2018,2021 Alessandro Sciarra
 * Copyright (c) 2014 Christopher Pinke
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

#include "benchmarkExecutable.hpp"

#include <numeric>

benchmarkExecutable::benchmarkExecutable(int argc, const char* argv[], std::string name)
    : generalExecutable(argc, argv, "benchmark"), kernelName(name)
{
    initializationTimer.reset();
    gaugefield = new physics::lattices::Gaugefield(*system,
                                                   &(interfacesHandler->getInterface<physics::lattices::Gaugefield>()),
                                                   *prng);
    /*
     * NOTE: Profiling is known to slow down the code in production, since a synchronization is required
     *       after kernel execution. It might be confusing to see this activated and used here. However,
     *       when running a benchmark on a single device, it is possible to take advantage of profiling to
     *       really measure execution time on device, getting rid of other host operations. For multiple
     *       device benchmarks, instead, this is not possible, since communication between devices is needed,
     *       and then we need to measure execution time on host. This is done using a timer and actually it is
     *       alsways done, so that the user get the host time information also when running on single device.
     *
     * NOTE: Changing the value of profiling here when the system has been already instantiated in the parent class
     *       can lead to inconsistencies, because the devices and in particular the OpenCL queue on them have to be
     *       initialized differently when profiling is on or off. We then fix the profiling variable depending on how
     *       many devices were requested/found in the previous system initialization and we then reinitialize system.
     *       It is worth remarking that hP and kP are aware of the change of enable_profiling since they contain a
     *       pointer to parameters as their member (parameter is a singleton in the codebase).
     */
    if (system->get_devices().size() == 1)
        parameters.enable_profiling = true;
    else
        parameters.enable_profiling = false;
    system         = new hardware::System(*hP, *kP);
    devices        = system->get_devices();
    benchmarkSteps = parameters.get_benchmarksteps();
    initializationTimer.add();
}

void benchmarkExecutable::runBenchmark()
{
    if (devices.size() > 1) {
        /*
         * Ensure kernels are built to avoid artificially deteriorating performance
         * This is not needed on single device, since performance is measured on device.
         */
        enqueueSpecificKernelForBenchmark();
    }
    logger.info() << "Performing benchmark...";
    performanceTimer.reset();
    for (unsigned int iteration = 0; iteration < benchmarkSteps; iteration++) {
        enqueueSpecificKernelForBenchmark();
    }
    performanceTimer.add();
    logger.info() << "...benchmark done!";
    printPerformanceToScreen();
}

void benchmarkExecutable::printPerformanceToScreen()
{
    auto kernelCallsNumber = 2 * benchmarkSteps;
    logger.info() << kernelName << " kernel called " << kernelCallsNumber << " times";
    double hostTime      = static_cast<double>(performanceTimer.getTime());
    double benchmarkTime = std::numeric_limits<double>::quiet_NaN();
    if (devices.size() == 1) {
        auto deviceTimes = getExecutionTimesOnDevices();
        benchmarkTime    = std::accumulate(deviceTimes.begin(), deviceTimes.end(), 0.0) / deviceTimes.size();
    } else {
        benchmarkTime = hostTime;
    }
    auto flopPerKernelCall   = getFlopsPerKernelCall();
    auto memoryPerKernelCall = getMemoryPerKernelCall();
    std::stringstream flopsToBePrinted, bytesToBePrinted;
    flopsToBePrinted << std::setprecision(3) << std::fixed << flopPerKernelCall / 1.e6 << " MFLOP";
    bytesToBePrinted << std::setprecision(3) << std::fixed << memoryPerKernelCall / 1.e6 << " MB";
    logger.info() << kernelName << " data per call: " << flopsToBePrinted.str() << " and " << bytesToBePrinted.str()
                  << " read/write memory";
    if (devices.size() == 1)
        logger.info() << kernelName << " device time per call: " << benchmarkTime / (kernelCallsNumber * 1.e3) << " ms";
    logger.info() << kernelName << " host time per call: " << hostTime / (kernelCallsNumber * 1.e3) << " ms";
    if (devices.size() == 1)
        logger.info() << "Benchmark done on single device, basing performance on device time:";
    else
        logger.info() << "Benchmark done on multiple devices, basing performance on host time:";
    double gflops = static_cast<double>(flopPerKernelCall) * kernelCallsNumber / benchmarkTime / 1e3;
    double gbytes = static_cast<double>(memoryPerKernelCall) * kernelCallsNumber / benchmarkTime / 1e3;
    logger.info() << kernelName << " performance: " << gflops << " GFLOPS";
    logger.info() << kernelName << " memory: " << gbytes << " GB/S";
}
