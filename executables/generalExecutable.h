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
 * Declaration of the generalExecutable class.
 * This class provides general features needed by the executables.
 */

#ifndef GENERALEXECUTABLE_H_
#define GENERALEXECUTABLE_H_

#include <fstream>

#include "../common_header_files/types.h"
#include "../host_functionality/host_use_timer.h"
#include "../executables/exceptions.h"
#include "../host_functionality/logger.hpp"
#include "../meta/util.hpp"
#include "../hardware/system.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/prng.hpp"
#include "../physics/observables/gaugeObservables.hpp"
#include "../physics/interfacesHandler.hpp"
#include "../interfaceImplementations/hardwareParameters.hpp"
#include "../interfaceImplementations/openClKernelParameters.hpp"

class generalExecutable
{

public:
	/**
	 * Destructor. Prints runtime information to screen and file.
	 */
	virtual ~generalExecutable();

protected:
	/**
	* Initialize meta::Inputparametes and hardware::System objects
	 */
	//Protected since it makes no sense to allow the user to instantiate this class
	generalExecutable(int argc, const char* argv[], std::string parameterSet = "all parameters"); 

	const char* ownName;
	std::string filenameForLogfile;
	std::string filenameForProfilingData;
	usetimer totalRuntimeOfExecutable;
	usetimer initializationTimer;
	usetimer performanceTimer;
	const meta::Inputparameters parameters;
	const physics::PrngParametersInterface * prngParameters;
	const hardware::HardwareParametersImplementation * hP;
	const hardware::code::OpenClKernelParametersImplementation * kP;
	hardware::System * system;
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
	std::unique_ptr<physics::InterfacesHandler> interfacesHandler;
	std::ofstream outputToFile;
	const char* generalTimeOutputFilename = "general_time_output";

	void printRuntimeInformationToScreenAndFile();

	void printGeneralTimesToScreen();

	void printGeneralTimesToFile();

	void printParametersToScreenAndFile();

	void printProfilingDataToFile();
};

#endif /* GENERALEXECUTABLE_H_ */
