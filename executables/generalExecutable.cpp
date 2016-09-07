/*
 * Copyright 2012, 2013, 2015 Lars Zeidlewicz, Christopher Pinke,
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

#include "generalExecutable.h"
#include "../interfaceImplementations/interfacesHandler.hpp"

void generalExecutable::printParametersToScreenAndFile()
{
    logger.info() << "## Starting executable: " << ownName;
    print_info_global(parameters);
    print_info_configs_io(parameters);
    print_info_prng_io(parameters);

    outputToFile.open(filenameForLogfile);
    if(outputToFile.is_open()) {
        meta::print_info_global(&outputToFile, parameters);
        meta::print_info_configs_io(&outputToFile, parameters);
        meta::print_info_prng_io(&outputToFile, parameters);
        outputToFile.close();
    } else {
        throw File_Exception(filenameForLogfile);
    }
}

void generalExecutable::printProfilingDataToFile()
{
    if(parameters.get_enable_profiling()) {
        logger.info() << "## writing general times to file: \"" << filenameForProfilingData << "\"";
        // For benchmarking one might need the lattice sizes
        outputToFile.open(filenameForProfilingData, std::ios::out | std::ios::app);
        if(outputToFile.is_open()) {
            meta::print_info_global(&outputToFile, parameters);
            outputToFile.close();
        } else {
            throw File_Exception(filenameForLogfile);
        }
        print_profiling(system, filenameForProfilingData);
    }
}

generalExecutable::generalExecutable(int argc, const char* argv[], std::string parameterSet)
        : parameters(argc, argv, parameterSet), prngParameters(nullptr)
{
	totalRuntimeOfExecutable.reset();
	initializationTimer.reset();
	ownName = argv[0];
	filenameForLogfile = meta::createLogfileName(ownName);
	filenameForProfilingData = meta::create_profiling_data_filename(parameters, ownName);
	switchLogLevel(parameters.get_log_level());
	printParametersToScreenAndFile();
	//@todo: these new here are not deleted apparently!!
	hP = new hardware::HardwareParametersImplementation (&parameters);
	kP = new hardware::code::OpenClKernelParametersImplementation (parameters);
	system = new hardware::System(*hP, *kP);
	prngParameters = new physics::PrngParametersImplementation(parameters);
	prng = new physics::PRNG(*system, prngParameters);
	interfacesHandler = std::unique_ptr<physics::InterfacesHandler>(new physics::InterfacesHandlerImplementation{parameters});
	initializationTimer.add();

}
generalExecutable::~generalExecutable()
{
    totalRuntimeOfExecutable.add();
    printRuntimeInformationToScreenAndFile();
    printProfilingDataToFile();
    if(prngParameters) {
        delete prngParameters;
    }
}

void generalExecutable::printRuntimeInformationToScreenAndFile()
{
    printGeneralTimesToScreen();
    printGeneralTimesToFile();
    return;
}

void generalExecutable::printGeneralTimesToScreen()
{
    using namespace std;
    logger.info() << "## *******************************************************************";
    logger.info() << "## General Times [mus]:";
    logger.info() << "## *******************************************************************";
    logger.info() << "## Program Parts:\t" << setfill(' ') << setw(5) << "total" << '\t' << setw(5) << "perc";
    logger.info() << "## Total:\t" << setfill(' ') << setw(12) << totalRuntimeOfExecutable.getTime();
    logger.info() << "## Init.:\t" << setfill(' ') << setw(12) << initializationTimer.getTime() << '\t' << fixed << setw(5) << setprecision(1)
            << percent(initializationTimer.getTime(), totalRuntimeOfExecutable.getTime());
    logger.info() << "## Perf.:\t" << setfill(' ') << setw(12) << performanceTimer.getTime() << '\t' << fixed << setw(5) << setprecision(1)
            << percent(performanceTimer.getTime(), totalRuntimeOfExecutable.getTime());
    logger.info() << "## *******************************************************************";
    return;
}
void generalExecutable::printGeneralTimesToFile()
{
    using namespace std;
    logger.info() << "## writing general times to file: \"" << generalTimeOutputFilename << "\"";
    outputToFile.open(generalTimeOutputFilename);
    if(outputToFile.is_open()) {
        outputToFile << "## *******************************************************************" << endl;
        outputToFile << "## General Times [mus]:" << endl;
        outputToFile << "## Total\tInit\tPerformance" << endl;
        outputToFile << totalRuntimeOfExecutable.getTime() << "\t" << initializationTimer.getTime() << '\t' << performanceTimer.getTime() << endl;
        outputToFile.close();
    } else {
        logger.warn() << "Could not open output file for general time output.";
    }
    return;
}

