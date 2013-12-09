/*
 * @file
 * Declaration of the hmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Hybrid Monte Carlo (HMC) algorithm.
 */

#ifndef HMCEXECUTABLE_H_
#define HMCEXECUTABLE_H_

#include "generationExecutable.h"
#include "physics/algorithms/hmc.hpp"
#include <cmath>

class hmcExecutable : public generationExecutable
{
public:
	hmcExecutable(int argc, const char* argv[]);

	~hmcExecutable();
protected:
	const std::string filenameForHmcLogfile = "hmc.log";
	double acceptanceRate = 0;
	hmc_observables observables;

	/*
	 * Sets member variables that control the iterations during
	 * the generation of gaugefield configurations.
	 */
	void setIterationParameters();

	void printParametersToScreenAndFile();

	void writeHmcLogfile();

	void thermalizeAccordingToSpecificAlgorithm();

	void generateAccordingToSpecificAlgorithm();

	void performOnlineMeasurements();

	void printHmcObservables(const std::string& filename);

	void printHmcObservablesToFile(const std::string& filename);

	void printHmcObservablesToScreen();
};

#endif /* HMCEXECUTABLE_H_ */
