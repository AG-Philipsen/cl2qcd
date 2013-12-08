/*
 * @file
 * Declaration of the hmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Hybrid Monte Carlo (HMC) algorithm.
 */

#ifndef HMCEXECUTABLE_H_
#define HMCEXECUTABLE_H_

#include "generationExecutable.h"

class hmcExecutable : public generationExecutable
{
public:
	hmcExecutable(int argc, const char* argv[]);

protected:
	/*
	 * Sets member variables that control the iterations during
	 * the generation of gaugefield configurations.
	 */
	void setIterationParameters();

	void thermalizeAccordingToSpecificAlgorithm();

	void generateAccordingToSpecificAlgorithm();

	void performOnlineMeasurements(int iteration);
};

#endif /* HMCEXECUTABLE_H_ */
