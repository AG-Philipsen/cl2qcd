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

#include "hmcExecutable.h"

#include "../physics/observables/wilsonTwoFlavourChiralCondensate.hpp"

hmcExecutable::hmcExecutable(int argc, const char* argv[])
        : generationExecutable(argc, argv, "hmc")
{
    initializationTimer.reset();
    printParametersToScreenAndFile();
    setIterationParameters();
    initializationTimer.add();
}

hmcExecutable::~hmcExecutable()
{
    using namespace std;
    logger.info() << "Acceptance rate: " << fixed << setprecision(1) << percent(acceptanceRate, parameters.get_hmcsteps()) << "%";
}

void hmcExecutable::printParametersToScreenAndFile()
{
    meta::print_info_hmc(parameters);
    writeHmcLogfile();
}

void hmcExecutable::writeHmcLogfile()
{
    outputToFile.open(filenameForLogfile, std::ios::out | std::ios::app);
    if(outputToFile.is_open()) {
        meta::print_info_hmc(&outputToFile, parameters);
        outputToFile.close();
    } else {
        throw File_Exception(filenameForLogfile);
    }
}

void hmcExecutable::setIterationParameters()
{
    generationExecutable::setIterationParameters();
    generationSteps += parameters.get_hmcsteps();
}

void hmcExecutable::thermalizeAccordingToSpecificAlgorithm()
{
    throw Print_Error_Message("Thermalization is not yet implemented for HMC algorithm");
}

void hmcExecutable::generateAccordingToSpecificAlgorithm()
{
    const double randomNumber = prng->get_double();
    observables = physics::algorithms::perform_hmc_step(gaugefield, iteration, randomNumber, *prng, *system, *interfacesHandler);
    acceptanceRate += observables.accept;
}

void hmcExecutable::performOnlineMeasurements()
{
    if(((iteration + 1) % writeFrequency) == 0) {
        std::string gaugeout_name = meta::get_hmc_obs_file_name(parameters, "");
        printHmcObservables(gaugeout_name);
        if(parameters.get_measure_pbp()) {
            physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(gaugefield, iteration, *interfacesHandler);
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
    const hmc_float exp_deltaH = std::exp(observables.deltaH);
    outputToFile.open(filename.c_str(), std::ios::out | std::ios::app);
    if(!outputToFile.is_open())
        throw File_Exception(filename);
    outputToFile << iteration << "\t";
    outputToFile.width(8);
    outputToFile.precision(15);
    outputToFile << observables.plaq << "\t" << observables.tplaq << "\t" << observables.splaq;
    outputToFile << "\t" << observables.poly.re << "\t" << observables.poly.im << "\t"
            << sqrt(observables.poly.re * observables.poly.re + observables.poly.im * observables.poly.im);
    outputToFile << "\t" << observables.deltaH << "\t" << exp_deltaH << "\t" << observables.prob << "\t" << observables.accept;
    //print number of iterations used in inversions with full and force precision
    /**
     * @todo: The counters should be implemented once the solver class is used!"
     * until then, only write "0"!
     */
    int iter0 = 0;
    int iter1 = 0;
    outputToFile << "\t" << iter0 << "\t" << iter1;
    if(parameters.get_use_mp()) {
        outputToFile << "\t" << iter0 << "\t" << iter1;
    }
    if(meta::get_use_rectangles(parameters)) {
        outputToFile << "\t" << observables.rectangles;
    }
    outputToFile << std::endl;
    outputToFile.close();
}

void hmcExecutable::printHmcObservablesToScreen()
{
    logger.info() << "\tHMC [OBS]:\t" << iteration << std::setw(8) << std::setfill(' ') << "\t" << std::setprecision(15) << observables.plaq << "\t"
            << observables.poly.re << "\t" << observables.poly.im;
}

