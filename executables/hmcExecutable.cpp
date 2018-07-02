/*
 * Copyright (c) 2013,2014 Christopher Pinke
 * Copyright (c) 2015,2018 Alessandro Sciarra
 * Copyright (c) 2015 Christopher Czaban
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

#include "hmcExecutable.hpp"

#include "../physics/observables/wilsonTwoFlavourChiralCondensate.hpp"

hmcExecutable::hmcExecutable(int argc, const char* argv[]) : generationExecutable(argc, argv, "hmc")
{
    initializationTimer.reset();
    printParametersToScreenAndFile();
    setIterationParameters();
    initializationTimer.add();
}

hmcExecutable::~hmcExecutable()
{
    using namespace std;
    logger.info() << "Acceptance rate: " << fixed << setprecision(1)
                  << percent(acceptanceRate, parameters.get_hmcsteps()) << "%";
}

void hmcExecutable::printParametersToScreenAndFile()
{
    meta::print_info_hmc(parameters);
    writeHmcLogfile();
}

void hmcExecutable::writeHmcLogfile()
{
    outputToFile.open(filenameForLogfile, std::ios::out | std::ios::app);
    if (outputToFile.is_open()) {
        meta::print_info_hmc(&outputToFile, parameters);
        outputToFile.close();
    } else {
        throw File_Exception(filenameForLogfile);
    }
}

void hmcExecutable::setIterationParameters()
{
    generationExecutable::setIterationParameters();
    nextToLastGenerationTraj += parameters.get_hmcsteps();
}

void hmcExecutable::thermalizeAccordingToSpecificAlgorithm()
{
    throw Print_Error_Message("Thermalization is not yet implemented for HMC algorithm");
}

void hmcExecutable::generateAccordingToSpecificAlgorithm()
{
    const double randomNumber = prng->get_double();
    observables = physics::algorithms::perform_hmc_step(gaugefield, iteration, randomNumber, *prng, *system,
                                                        *interfacesHandler);
    acceptanceRate += observables.accept;
}

void hmcExecutable::performOnlineMeasurements()
{
    if (((iteration + 1) % writeFrequency) == 0) {
        std::string gaugeout_name = meta::get_hmc_obs_file_name(parameters, "");
        printHmcObservables(gaugeout_name);
        if (parameters.get_measure_pbp()) {
            physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(gaugefield, iteration,
                                                                                          *interfacesHandler);
        }
    }
}

void hmcExecutable::printHmcObservables(const std::string& filename)
{
    printHmcObservablesToFile(filename);
    printHmcObservablesToScreen();
}

void hmcExecutable::printHmcObservablesToFile(const std::string& filename)
{
    outputToFile.open(filename.c_str(), std::ios::out | std::ios::app);
    if (!outputToFile.is_open())
        throw File_Exception(filename);
    const std::streamsize shortPrecision = 4;
    const std::streamsize longPrecision  = 15;
    const std::streamsize shortWidth     = shortPrecision + 6;
    const std::streamsize longWidth      = longPrecision + 10;  // +1 is always needed for the period, +4 is for e+XX
                                                                // in case of extreme values, +10 to give some breath
    outputToFile.precision(longPrecision);
    outputToFile << std::setw(8) << iteration  // statistics up to 1e8-1
                 << ' ' << std::setw(longWidth) << observables.plaq << ' ' << std::setw(longWidth) << observables.tplaq
                 << ' ' << std::setw(longWidth) << observables.splaq << ' ' << std::setw(longWidth)
                 << observables.poly.re << ' ' << std::setw(longWidth) << observables.poly.im << ' '
                 << std::setw(longWidth)
                 << sqrt(observables.poly.re * observables.poly.re + observables.poly.im * observables.poly.im) << ' '
                 << std::setw(longWidth) << observables.deltaH;
    outputToFile.precision(shortPrecision);
    outputToFile << ' ' << std::setw(6)
                 << observables.accept  // we print 0 or 1 with some space around, but not too much
                 << ' ' << std::setw(shortWidth) << observables.timeTrajectory;

    /**
     * @TODO: Add here to the files the number of iterations used in inversions with high and low precision.
     *        The counters should be implemented once the solver class is used! Something like:
     *            int iter0 = 0;
     *            int iter1 = 0;
     *            outputToFile << "\t" << iter0 << "\t" << iter1;
     *            if(parameters.get_use_mp()) {
     *                outputToFile << "\t" << iter0 << "\t" << iter1;
     *            }
     */
    if (meta::get_use_rectangles(parameters)) {
        outputToFile.precision(longPrecision);
        outputToFile << ' ' << std::setw(longWidth) << observables.rectangles;
    }
    outputToFile << std::endl;
    outputToFile.close();
}

void hmcExecutable::printHmcObservablesToScreen()
{
    logger.info() << "\tHMC [OBS]:\t" << iteration << std::setw(8) << std::setfill(' ') << "\t" << std::setprecision(15)
                  << observables.plaq << "\t" << observables.poly.re << "\t" << observables.poly.im;
}
