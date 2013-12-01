/*
 * @file
 * Declaration of the generalExecutable class.
 * This class provides general features needed by the executables.
 */

#ifndef GENERALEXECUTABLE_H_
#define GENERALEXECUTABLE_H_

#include <fstream>

#include "types.h"
#include "host_use_timer.h"
#include "exceptions.h"
#include "logger.hpp"
#include "meta/util.hpp"
#include "./hardware/system.hpp"

class generalExecutable
{

public:
	/**
	 * Initialize meta::Inputparametes and Shardware::System objects
	 */
	generalExecutable(int argc, const char* argv[]);

	/**
	 * Destructor. Prints runtime information to screen and file.
	 */
	~generalExecutable();

protected:
	const char* ownName;
	usetimer totalRuntimeOfExecutable;
	usetimer initializationTimer;
	usetimer performanceTimer;
	meta::Inputparameters parameters;
	hardware::System * system;
	std::ofstream outputToFile;
	const char* generalTimeOutputFilename = "general_time_output";

	void printRuntimeInformationToScreenAndFile();

	void printGeneralTimesToScreen();

	void printGeneralTimesToFile();

};

#endif /* GENERALEXECUTABLE_H_ */
