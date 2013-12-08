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

protected:
	const std::string filenameForHmcLogfile = "hmc.log";
	double acceptanceRate = 0;
	hmc_observables observables;

	/*
	 * Sets member variables that control the iterations during
	 * the generation of gaugefield configurations.
	 */
	void setIterationParameters();

	void writeHmcLogfile();

	void thermalizeAccordingToSpecificAlgorithm();

	void generateAccordingToSpecificAlgorithm();

	void performOnlineMeasurements(int iteration);

	void print_hmcobservables(const hmc_observables& obs, int iter, const std::string& filename, const meta::Inputparameters& params);

	void print_hmcobservables(const hmc_observables& obs, int iter);
};

#endif /* HMCEXECUTABLE_H_ */
